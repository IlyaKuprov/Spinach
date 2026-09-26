# examples/benchmarks/fft_gpu.m

- Signature: `fft_gpu()`

## Purpose

GPU arithmetic benchmark -3D Fourier transforms.

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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
