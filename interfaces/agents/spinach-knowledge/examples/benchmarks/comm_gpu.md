# examples/benchmarks/comm_gpu.m

- Signature: `comm_gpu(n)`

## Purpose

GPU communications benchmark. Adapted from example code in Matlab documentation. Data for IK's favourite NVIDIA cards on Dell Power Edge T630 workstation: Tesla K40 (2013): send 9.7, gather 2.6, bw GPU 190, bw CPU 64 Titan V (2017, TCC mode): send 10.3, gather 2.6, bw GPU 568, bw CPU 64 Tesla A100 (2021, PCI-E): send 10.4, gather 2.6, bw GPU 1291, bw CPU 64

## Implementation structure

- GPU communications benchmark. Adapted from example code in Matlab
- documentation. Data for IK's favourite NVIDIA cards on Dell Power
- Edge T630 workstation:
- Tesla K40 (2013): send 9.7, gather 2.6, bw GPU 190, bw CPU 64
- Titan V (2017, TCC mode): send 10.3, gather 2.6, bw GPU 568, bw CPU 64
- Tesla A100 (2021, PCI-E): send 10.4, gather 2.6, bw GPU 1291, bw CPU 64
- All GPUs by default
- Pick the GPU
- 8 bytes per double
- Array sizes to test
- Preallocate answer arrays
- Measure performance
