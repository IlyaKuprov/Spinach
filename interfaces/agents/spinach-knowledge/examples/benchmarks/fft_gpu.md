# examples/benchmarks/fft_gpu.m

- MATLAB implementation: [examples/benchmarks/fft_gpu.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/benchmarks/fft_gpu.m)

- Signature: `fft_gpu()`.
- Returns: no MATLAB output arguments. Prints per-run CPU/GPU timings and opens a comparison plot.

## Purpose and use

This script benchmarks MATLAB's three-dimensional `fftn` on cubic double-precision random arrays, comparing a CPU array with a `gpuArray` on the default GPU device. Call `fft_gpu()`; it takes no configuration arguments. If `gpuDeviceCount` reports zero, it prints `no CUDA GPUs detected` and returns without benchmarking. Otherwise it initialises `gpuDevice` and uses that device. MATLAB GPU support and enough device memory for the chosen arrays are required.

## Benchmark procedure

The cube side lengths are 128, 192, 256, 384, and 512 points. For each size, ten CPU and ten GPU runs are timed with `tic`/`toc`; each iteration creates a new `randn(...,'double')` CPU array and a `randn(...,'gpuArray')` GPU array before timing its `fftn`. `wait(dev)` brackets the GPU transform so asynchronous device work is complete at both timing boundaries. The script prints every individual timing in seconds.

The plot uses the mean of repetitions 2–10 (the first repetition is excluded) for each implementation. It plots CPU and GPU means against the cube side length—not total voxel count—with logarithmic axes, labels the ordinate as calculation time in seconds, and restricts the displayed x range to 100–600. The first sample is discarded from the displayed mean, but the source does not explain why.

## Limits

This measures FFT throughput, not a Spinach spin-dynamics workload or an end-to-end scientific calculation. The largest cubes can exhaust GPU memory; the source comment recommends reducing `sizes` if the card runs out. Timings and the speed comparison are hardware- and MATLAB-dependent, and the script does not return the timing array or save the figure automatically.
