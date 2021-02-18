/*
 * sketchlib_bindings.cpp
 * Python bindings for pp-sketchlib
 *
 */

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

constexpr cuda_version = __CUDACC_VER_MAJOR__ + "."
                         __CUDACC_VER_MINOR__ + "."
                         __CUDACC_VER_BUILD__;

PYBIND11_MODULE(sparrowhawk, m)
{
  m.doc() = "Microbial genome assembly using CUDA";

  // Exported functions
  m.def("deviceInfo", &deviceInfo, "Get information about available GPUs");

  m.def("makeDBG", &makeDBG, "Make a DBG at a given k-mer length fro sequence reads",
        py::arg("fwd"),
        py::arg("rev"),
        py::arg("k"),
        py::arg("output"),
        py::arg("num_threads") = 1,
        py::arg("device_id") = 0);

  m.def("assembleReads", &assembleReads, "Assemble contigs from sequence reads",
        py::arg("fwd"),
        py::arg("rev"),
        py::arg("klist"),
        py::arg("output"),
        py::arg("num_threads") = 1,
        py::arg("device_id") = 0);

  m.attr("version") = VERSION_INFO;
  m.attr("cudaVersion") = cuda_version;
}
