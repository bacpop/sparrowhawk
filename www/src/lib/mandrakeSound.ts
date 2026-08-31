/**
 * Live sonification of the Mandrake embedding as it converges. Modelled on
 * mandrake/sound.py + src/sound.hpp from the original C++/Python tool, which
 * baked a similar tone sequence into an offline --animate video. Here it plays
 * directly through the Web Audio API as "frame" updates arrive from the
 * worker, instead of synthesizing one full sample buffer up front.
 */

const FREQ_MIN = 120;
const FREQ_MAX = 1200;
const TONE_DURATION = 0.18; // seconds
const ATTACK = 0.025;
const DECAY = 0.1;
const SUSTAIN = 0.9;
const RELEASE = 0.3;
const PEAK_GAIN = 0.12;

interface Frame {
  x: Float64Array;
  y: Float64Array;
}

// Coordinates drift in absolute scale/position during optimisation (the
// whole layout can grow, shrink, or shift, especially with initial
// exaggeration), so frames are centred and rescaled before comparing them —
// otherwise a global scale change would swamp genuine point movement.
function centreAndScale(embedding: Float64Array): Frame {
  const pointCount = Math.floor(embedding.length / 2);
  const x = new Float64Array(pointCount);
  const y = new Float64Array(pointCount);
  let meanX = 0;
  let meanY = 0;
  for (let index = 0; index < pointCount; index += 1) {
    meanX += embedding[index * 2];
    meanY += embedding[index * 2 + 1];
  }
  meanX /= pointCount || 1;
  meanY /= pointCount || 1;
  let maxAbs = 0;
  for (let index = 0; index < pointCount; index += 1) {
    x[index] = embedding[index * 2] - meanX;
    y[index] = embedding[index * 2 + 1] - meanY;
    maxAbs = Math.max(maxAbs, Math.abs(x[index]), Math.abs(y[index]));
  }
  if (maxAbs > 0) {
    for (let index = 0; index < pointCount; index += 1) {
      x[index] /= maxAbs;
      y[index] /= maxAbs;
    }
  }
  return { x, y };
}

// Per-axis maximum absolute displacement between two frames — the same
// convergence measure sound.py computes offline (np.max(np.abs(delta),
// axis=0)) before mapping it to pitch. A still layout stays silent-ish
// (low, unchanging pitch); a layout still rearranging itself sounds active.
function maxAbsDelta(previous: Float64Array, next: Float64Array): number {
  let max = 0;
  const count = Math.min(previous.length, next.length);
  for (let index = 0; index < count; index += 1) {
    max = Math.max(max, Math.abs(next[index] - previous[index]));
  }
  return max;
}

/** Envelope shape ported from src/sound.hpp's oscillator::envelope. */
function scheduleEnvelope(gain: GainNode, startTime: number, duration: number, peak: number): void {
  const attackEnd = startTime + duration * ATTACK;
  const decayEnd = attackEnd + duration * DECAY;
  const releaseStart = startTime + duration * (1 - RELEASE);
  const end = startTime + duration;
  gain.gain.cancelScheduledValues(startTime);
  gain.gain.setValueAtTime(0, startTime);
  gain.gain.linearRampToValueAtTime(peak, attackEnd);
  gain.gain.linearRampToValueAtTime(peak * SUSTAIN, decayEnd);
  gain.gain.setValueAtTime(peak * SUSTAIN, releaseStart);
  gain.gain.linearRampToValueAtTime(0, end);
}

function playTone(context: AudioContext, frequency: number, pan: number): void {
  const oscillator = context.createOscillator();
  oscillator.type = "triangle";
  oscillator.frequency.value = frequency;
  const gain = context.createGain();
  const panner = context.createStereoPanner();
  panner.pan.value = pan;
  oscillator.connect(gain).connect(panner).connect(context.destination);
  const startTime = context.currentTime;
  scheduleEnvelope(gain, startTime, TONE_DURATION, PEAK_GAIN);
  oscillator.start(startTime);
  oscillator.stop(startTime + TONE_DURATION + 0.02);
  oscillator.onended = () => {
    oscillator.disconnect();
    gain.disconnect();
    panner.disconnect();
  };
}

// sound.py knows every frame's displacement up front (the run already
// finished) and normalises against the batch min/max. Live frames arrive
// one at a time, so this tracks an incrementally-widening min/max instead —
// pitch mapping stabilises as more frames are seen within the same run.
class RunningRange {
  private min = Infinity;
  private max = 0;

  toFrequency(value: number): number {
    this.min = Math.min(this.min, value);
    this.max = Math.max(this.max, value);
    const span = this.max - this.min;
    const fraction = span > 0 ? (value - this.min) / span : 0.5;
    return FREQ_MIN + fraction * (FREQ_MAX - FREQ_MIN);
  }

  reset(): void {
    this.min = Infinity;
    this.max = 0;
  }
}

export class MandrakeSonifier {
  private context: AudioContext | null = null;
  private previous: Frame | null = null;
  private xRange = new RunningRange();
  private yRange = new RunningRange();

  /** Create/resume the AudioContext. Call from a user gesture. */
  enable(): void {
    if (!this.context) this.context = new AudioContext();
    if (this.context.state === "suspended") void this.context.resume();
    this.previous = null;
    this.xRange.reset();
    this.yRange.reset();
  }

  disable(): void {
    this.previous = null;
  }

  dispose(): void {
    void this.context?.close();
    this.context = null;
    this.previous = null;
  }

  /** Play a tone pair for the movement between the previous and this frame. */
  playFrame(embedding: Float64Array): void {
    if (!this.context) return;
    const current = centreAndScale(embedding);
    if (this.previous) {
      const dx = maxAbsDelta(this.previous.x, current.x);
      const dy = maxAbsDelta(this.previous.y, current.y);
      playTone(this.context, this.xRange.toFrequency(dx), -0.6);
      playTone(this.context, this.yRange.toFrequency(dy), 0.6);
    }
    this.previous = current;
  }
}
