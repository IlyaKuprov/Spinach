# examples/benchmarks/comm_gpu.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/comm_gpu.m`
- Signature: `comm_gpu(n)`
- Total lines: 97

## Purpose

GPU communications benchmark. Adapted from example code in Matlab documentation. Data for IK's favourite NVIDIA cards on Dell Power Edge T630 workstation: Tesla K40 (2013): send 9.7, gather 2.6, bw GPU 190, bw CPU 64 Titan V (2017, TCC mode): send 10.3, gather 2.6, bw GPU 568, bw CPU 64 Tesla A100 (2021, PCI-E): send 10.4, gather 2.6, bw GPU 1291, bw CPU 64

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: All GPUs by default; implemented by `if nargin==0`.
- Lines 21-22: Pick the GPU; implemented by `gpu=gpuDevice(n)`.
- Lines 26-27: 8 bytes per double; implemented by `size_of_double=8`.
- Lines 29-30: Array sizes to test; implemented by `sizes=power(2,14:28)`.
- Lines 32-33: Preallocate answer arrays; implemented by `send_times=inf(size(sizes))`.
- Lines 38-39: Measure performance; implemented by `for n=1:numel(sizes)`.
- Lines 41-42: Generate arrays; implemented by `num_elements=sizes(n)/size_of_double`.
- Lines 46-47: Time sending to GPU; implemented by `send_fcn=@()gpuArray(host_data)`.
- Lines 50-51: Time gathering from GPU; implemented by `gather_fcn=@()gather(gpu_data)`.
- Lines 54-55: Time read+write inside GPU; implemented by `plus_fcn_gpu=@()plus(gpu_data,1.0)`.
- Lines 58-59: Time read+write on host; implemented by `plus_fcn_host=@()plus(host_data,1.0)`.
- Lines 64-65: Get performance figures; implemented by `send_bandwidth=(sizes./send_times)/1e9`.
- Lines 74-75: Report to console; implemented by `disp(['Peak send bandwidth: ' num2str(max_send_bandwidth) ' GB/s'])`.
- Lines 80-81: Do the plotting; implemented by `figure('Name',namestring,'NumberTitle','off')`.

### Control flow inferred from the code

- Line 14: conditional branch on `nargin==0`.
- Line 15: `for` loop over `n=1:gpuDeviceCount('available')`.
- Line 39: `for` loop over `n=1:numel(sizes)`.

### Key state/data transformations

- Lines 22: computes `gpu` using `gpu=gpuDevice(n)`.
- Lines 23: computes `namestring` using `namestring=['GPU ' num2str(n) ': ' gpu.Name]`.
- Lines 27: computes `size_of_double` using `size_of_double=8`.
- Lines 30: computes `sizes` using `sizes=power(2,14:28)`.
- Lines 33: computes `send_times` using `send_times=inf(size(sizes))`.
- Lines 34: computes `gather_times` using `gather_times=inf(size(sizes))`.
- Lines 35: computes `memory_times_gpu` using `memory_times_gpu=inf(size(sizes))`.
- Lines 36: computes `memory_times_host` using `memory_times_host=inf(size(sizes))`.
- Lines 42: computes `num_elements` using `num_elements=sizes(n)/size_of_double`.
- Lines 43: computes `host_data` using `host_data=randi([0 9],num_elements,1)`.
- Lines 44: computes `gpu_data` using `gpu_data=randi([0 9],num_elements,1,'gpuArray')`.
- Lines 47: computes `send_fcn` using `send_fcn=@()gpuArray(host_data)`.
- Lines 48: computes `send_times(n)` using `send_times(n)=gputimeit(send_fcn)`.
- Lines 51: computes `gather_fcn` using `gather_fcn=@()gather(gpu_data)`.
- Lines 52: computes `gather_times(n)` using `gather_times(n)=gputimeit(gather_fcn)`.
- Lines 55: computes `plus_fcn_gpu` using `plus_fcn_gpu=@()plus(gpu_data,1.0)`.
- Lines 56: computes `memory_times_gpu(n)` using `memory_times_gpu(n)=gputimeit(plus_fcn_gpu)`.
- Lines 59: computes `plus_fcn_host` using `plus_fcn_host=@()plus(host_data,1.0)`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gpuDeviceCount()`, `gpuDevice()`, `num2str()`, `power()`, `inf()`, `sizes()`, `randi()`, `gpuArray()`, `send_times()`, `gputimeit()`, `gather()`, `gather_times()`, `plus()`, `memory_times_gpu()`, `memory_times_host()`, `timeit()`.
