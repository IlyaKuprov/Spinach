# examples/benchmarks/mult_gpu.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/mult_gpu.m`
- Signature: `mult_gpu(precision)`
- Total lines: 69

## Purpose

CPU and GPU matrix arithmetic benchmark. Set the argument to either 'single' or 'double'. IK's workstation output: GPU, precision: single Titan V PCIe (2017): 12 TFLOPS Tesla A100 PCIe (2021): 10 TFLOPS Tesla A800 PCIe (2024): 18 TFLOPS Tesla H200 SXM (2025): 51 TFLOPS GPU, precision: double Titan V PCIe (2017): 6.5 TFLOPS Tesla A100 PCIe (2021): 10 TFLOPS Tesla A800 PCIe (2024): 15 TFLOPS Tesla H200 SXM (2025): 60 T

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Default to double precision; implemented by `if ~exist('precision','var')`.
- Lines 26-27: Look for GPUs; implemented by `if gpuDeviceCount==0`.
- Lines 31-32: Set matrix sizes; implemented by `sizes=power(2,20:2:26); N=sqrt(sizes)`.
- Lines 34-35: Preallocate timing arrays; implemented by `mmTimesCPU=nan(size(sizes))`.
- Lines 38-39: Loop over sizes; implemented by `for n=1:numel(sizes)`.
- Lines 41-42: Run CPU benchmark; implemented by `A=rand(N(n),N(n),precision)`.
- Lines 46-47: Run GPU benchmark; implemented by `if gpuDeviceCount>0`.
- Lines 54-55: Convert to GFlops; implemented by `mmTFlopsCPU=(2*N.^3-N.^2)./mmTimesCPU/1e12`.
- Lines 58-59: Get the maxima; implemented by `maxTFlopsCPU=max(mmTFlopsCPU)`.
- Lines 62-63: Report to the user; implemented by `disp(['Precision: ' precision])`.

### Control flow inferred from the code

- Line 22: conditional branch on `~exist('precision','var')`.
- Line 27: conditional branch on `gpuDeviceCount==0`.
- Line 39: `for` loop over `n=1:numel(sizes)`.
- Line 47: conditional branch on `gpuDeviceCount>0`.

### Key state/data transformations

- Lines 23: computes `precision` using `precision='double'`.
- Lines 32: computes `sizes` using `sizes=power(2,20:2:26); N=sqrt(sizes)`.
- Lines 35: computes `mmTimesCPU` using `mmTimesCPU=nan(size(sizes))`.
- Lines 36: computes `mmTimesGPU` using `mmTimesGPU=nan(size(sizes))`.
- Lines 42: computes `A` using `A=rand(N(n),N(n),precision)`.
- Lines 43: computes `B` using `B=rand(N(n),N(n),precision)`.
- Lines 44: computes `mmTimesCPU(n)` using `mmTimesCPU(n)=timeit(@()A*B)`.
- Lines 49: computes `mmTimesGPU(n)` using `mmTimesGPU(n)=gputimeit(@()A*B)`.
- Lines 55: computes `mmTFlopsCPU` using `mmTFlopsCPU=(2*N.^3-N.^2)./mmTimesCPU/1e12`.
- Lines 56: computes `mmTFlopsGPU` using `mmTFlopsGPU=(2*N.^3-N.^2)./mmTimesGPU/1e12`.
- Lines 59: computes `maxTFlopsCPU` using `maxTFlopsCPU=max(mmTFlopsCPU)`.
- Lines 60: computes `maxTFlopsGPU` using `maxTFlopsGPU=max(mmTFlopsGPU)`.

## Implementation structure

- CPU and GPU matrix arithmetic benchmark. Set the argument to
- either 'single' or 'double'. IK's workstation output:
- GPU, precision: single
- Titan V PCIe (2017): 12 TFLOPS
- Tesla A100 PCIe (2021): 10 TFLOPS
- Tesla A800 PCIe (2024): 18 TFLOPS
- Tesla H200 SXM (2025): 51 TFLOPS
- GPU, precision: double
- Titan V PCIe (2017): 6.5 TFLOPS
- Tesla A800 PCIe (2024): 15 TFLOPS
- Tesla H200 SXM (2025): 60 TFLOPS
- Default to double precision

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `power()`, `nan()`, `mmTimesCPU()`, `timeit()`, `gpuArray()`, `mmTimesGPU()`, `gputimeit()`, `TFLOPS()`.
