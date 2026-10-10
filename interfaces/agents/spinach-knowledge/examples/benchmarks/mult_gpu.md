# examples/benchmarks/mult_gpu.m

## Use

mult_gpu(precision) compares dense square matrix-multiplication throughput on the CPU and, when a CUDA GPU is detected, the GPU. Call mult_gpu() or mult_gpu('double') for the default double precision, or mult_gpu('single'). The source documents single and double; it does not validate the argument itself. The function returns no values and prints the precision and best reported CPU/GPU rates.

## Benchmark method and limits

It benchmarks random square matrices with dimensions 1,024, 2,048, 4,096, and 8,192. CPU timings use timeit; if gpuDeviceCount is positive, the matrices are copied to GPU arrays and timed with gputimeit. For each size it converts elapsed time to TFLOPS using the operation count 2*N^3-N^2 divided by seconds and 1e12, then reports the maximum across sizes. The largest case uses two 8,192-by-8,192 matrices, so sufficient host and (for GPU timing) device memory is needed.

With no detected GPU, it displays “no CUDA GPUs detected” and skips GPU timing; the GPU timing entries remain initialised to NaN. This is a benchmark, not a device-independent performance guarantee. The hardware figures below are the workstation results recorded in the source comments, not expected results for other machines.

## Recorded workstation results

| GPU | Precision | Recorded throughput |
|---|---|---:|
| Titan V PCIe (2017) | single | 12 TFLOPS |
| Tesla A100 PCIe (2021) | single | 10 TFLOPS |
| Tesla A800 PCIe (2024) | single | 18 TFLOPS |
| Tesla H200 SXM (2025) | single | 51 TFLOPS |
| Titan V PCIe (2017) | double | 6.5 TFLOPS |
| Tesla A100 PCIe (2021) | double | 10 TFLOPS |
| Tesla A800 PCIe (2024) | double | 15 TFLOPS |
| Tesla H200 SXM (2025) | double | 60 TFLOPS |

## Source

[examples/benchmarks/mult_gpu.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/benchmarks/mult_gpu.m)
