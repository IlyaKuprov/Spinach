# examples/benchmarks/fft_gpu.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/fft_gpu.m`
- Signature: `fft_gpu()`
- Total lines: 54

## Purpose

GPU arithmetic benchmark -3D Fourier transforms.

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Look for GPUs; implemented by `if gpuDeviceCount==0`.
- Lines 12-13: Initialise GPU device; implemented by `dev=gpuDevice`.
- Lines 15-16: FFT dimensions (reduce if card runs out of memory); implemented by `sizes=[128 192 256 384 512]`.
- Lines 18-19: Timing array; implemented by `timings=zeros(numel(sizes),2,10)`.
- Lines 21-22: FFT size loop; implemented by `for n=1:numel(sizes)`.
- Lines 24-25: Statistics loop; implemented by `for k=1:10`.
- Lines 27-28: CPU benchmark; implemented by `a=randn(sizes(n),sizes(n),sizes(n),'double')`.
- Lines 32-33: GPU benchmark; implemented by `b=randn(sizes(n),sizes(n),sizes(n),'gpuArray')`.
- Lines 41-42: Analysis; implemented by `means=mean(timings(:,:,2:end),3)`.
- Lines 44-45: Plotting; implemented by `kfigure(); hold on`.

### Control flow inferred from the code

- Line 8: conditional branch on `gpuDeviceCount==0`.
- Line 22: `for` loop over `n=1:numel(sizes)`.
- Line 25: `for` loop over `k=1:10`.

### Key state/data transformations

- Lines 13: computes `dev` using `dev=gpuDevice`.
- Lines 16: computes `sizes` using `sizes=[128 192 256 384 512]`.
- Lines 19: computes `timings` using `timings=zeros(numel(sizes),2,10)`.
- Lines 28: computes `a` using `a=randn(sizes(n),sizes(n),sizes(n),'double')`.
- Lines 29: computes `tic; fftn(a); timings(n,1,k)` using `tic; fftn(a); timings(n,1,k)=toc`.
- Lines 30: computes `disp(['n` using `disp(['n=' num2str(sizes(n)) ', CPU time ' num2str(timings(n,1,k)) ' seconds'])`.
- Lines 33: computes `b` using `b=randn(sizes(n),sizes(n),sizes(n),'gpuArray')`.
- Lines 34: computes `wait(dev); tic; fftn(b); wait(dev); timings(n,2,k)` using `wait(dev); tic; fftn(b); wait(dev); timings(n,2,k)=toc`.
- Lines 42: computes `means` using `means=mean(timings(:,:,2:end),3)`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sizes()`, `timings()`, `num2str()`, `wait()`, `kfigure()`, `means()`, `xlim()`, `set()`.
