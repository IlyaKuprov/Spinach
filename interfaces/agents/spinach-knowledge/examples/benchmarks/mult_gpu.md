# examples/benchmarks/mult_gpu.m

- Signature: `mult_gpu(precision)`

## Purpose

CPU and GPU matrix arithmetic benchmark. Set the argument to either 'single' or 'double'. IK's workstation output: GPU, precision: single Titan V PCIe (2017): 12 TFLOPS Tesla A100 PCIe (2021): 10 TFLOPS Tesla A800 PCIe (2024): 18 TFLOPS Tesla H200 SXM (2025): 51 TFLOPS GPU, precision: double Titan V PCIe (2017): 6.5 TFLOPS Tesla A100 PCIe (2021): 10 TFLOPS Tesla A800 PCIe (2024): 15 TFLOPS Tesla H200 SXM (2025): 60 T

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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
