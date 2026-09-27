# examples/benchmarks/fft_gpu.m

- Signature: `fft_gpu()`

## Purpose

GPU arithmetic benchmark -3D Fourier transforms.

## Implementation structure

- GPU arithmetic benchmark -3D Fourier transforms.
- Look for GPUs
- Initialise GPU device
- FFT dimensions (reduce if card runs out of memory)
- Timing array
- FFT size loop
- Statistics loop
- CPU benchmark
- GPU benchmark
- Analysis
- Plotting
