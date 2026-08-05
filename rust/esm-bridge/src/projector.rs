//! The 2-D UMAP encoder, trained by `src/bin/train_umap.rs` and applied here.

use burn::nn::{LeakyRelu, LeakyReluConfig, Linear, LinearConfig};
use burn::prelude::*;
use burn::tensor::Bytes;
use burn_store::{BurnpackStore, ModuleSnapshot};

use crate::HIDDEN;

/// Hidden layer widths. Must match the trained artefact, so the trainer reads it from here.
pub const HIDDEN_LAYERS: [usize; 2] = [256, 128];
pub const N_DIM: usize = 2;

const NEGATIVE_SLOPE: f64 = 0.1;

/// The encoder. One definition, used by `train_umap fit` and by the browser alike, so there
/// is no architecture to keep in step between them.
#[derive(Module, Debug)]
pub struct UmapMlp<B: Backend> {
    layers: Vec<Linear<B>>,
    activation: LeakyRelu,
}

impl<B: Backend> UmapMlp<B> {
    pub fn new(device: &B::Device) -> Self {
        let mut layers = Vec::with_capacity(HIDDEN_LAYERS.len() + 1);
        let mut input = HIDDEN;
        for &width in &HIDDEN_LAYERS {
            layers.push(LinearConfig::new(input, width).with_bias(true).init(device));
            input = width;
        }
        layers.push(LinearConfig::new(input, N_DIM).with_bias(true).init(device));
        Self {
            layers,
            activation: LeakyReluConfig::new().with_negative_slope(NEGATIVE_SLOPE).init(),
        }
    }

    /// Load a burnpack written by `train_umap`. `Err` means the tab runs without a plot.
    pub fn from_bytes(bytes: Vec<u8>, device: &B::Device) -> Result<Self, String> {
        let mut model = Self::new(device);
        let mut store = BurnpackStore::from_bytes(Some(Bytes::from_bytes_vec(bytes)));
        model
            .load_from(&mut store)
            .map_err(|e| format!("could not load the UMAP encoder: {e}"))?;
        Ok(model)
    }

    pub fn forward(&self, input: Tensor<B, 2>) -> Tensor<B, 2> {
        let mut x = input;
        for layer in &self.layers[..self.layers.len() - 1] {
            x = self.activation.forward(layer.forward(x));
        }
        // No activation on the output.
        self.layers.last().unwrap().forward(x)
    }

    /// Row-major `[n, HIDDEN]` in, row-major `[n, N_DIM]` out.
    pub fn project(&self, vectors: &[f32], device: &B::Device) -> Vec<f32> {
        let n = vectors.len() / HIDDEN;
        if n == 0 {
            return Vec::new();
        }
        let input = Tensor::<B, 2>::from_data(
            TensorData::new(vectors.to_vec(), [n, HIDDEN]),
            device,
        );
        self.forward(input)
            .into_data()
            .to_vec()
            .expect("projection output should be f32")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CpuDevice;

    #[test]
    fn projects_to_two_dimensions() {
        let device = CpuDevice::default();
        let model = UmapMlp::<crate::Cpu>::new(&device);
        let coords = model.project(&vec![0.5f32; HIDDEN * 3], &device);
        assert_eq!(coords.len(), 3 * N_DIM);
        assert!(coords.iter().all(|c| c.is_finite()));
    }

    #[test]
    fn an_empty_matrix_projects_to_nothing() {
        let device = CpuDevice::default();
        let model = UmapMlp::<crate::Cpu>::new(&device);
        assert!(model.project(&[], &device).is_empty());
    }

    /// Save and reload through the same path `train_umap fit` and the browser use.
    #[test]
    fn round_trips_through_a_burnpack() {
        use burn_store::ModuleSnapshot;
        let device = CpuDevice::default();
        let model = UmapMlp::<crate::Cpu>::new(&device);
        let input: Vec<f32> = (0..HIDDEN * 4).map(|i| (i as f32 * 0.01).sin()).collect();
        let before = model.project(&input, &device);

        let mut out = BurnpackStore::from_bytes(None);
        model.save_into(&mut out).expect("save");
        let bytes = out.get_bytes().expect("bytes").to_vec();

        let reloaded = UmapMlp::<crate::Cpu>::from_bytes(bytes, &device).expect("reload");
        assert_eq!(before, reloaded.project(&input, &device));
    }

    #[test]
    fn rejects_bytes_that_are_not_an_encoder() {
        let device = CpuDevice::default();
        assert!(UmapMlp::<crate::Cpu>::from_bytes(vec![0u8; 32], &device).is_err());
    }
}
