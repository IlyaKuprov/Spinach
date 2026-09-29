# examples/benchmarks/comm_gpu.m

- MATLAB implementation: [examples/benchmarks/comm_gpu.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/benchmarks/comm_gpu.m)

- Signature: `comm_gpu(n)`; with no argument, runs once for each GPU reported by `gpuDeviceCount('available')`.
- Returns: no MATLAB output arguments. Prints peak bandwidths and opens a two-panel figure.

## Purpose and use

This standalone MATLAB GPU benchmark separates host-to-device transfer, device-to-host transfer, and simple read/write throughput on the GPU and CPU. Run `comm_gpu()` to use every available device in turn, or `comm_gpu(n)` to select device index `n` through `gpuDevice(n)`. It requires MATLAB GPU support and an available CUDA GPU; it is not a Spinach simulation.

## What it measures

The test sizes are `2.^(14:28)` bytes. Each size is divided by eight to allocate a column of double-precision values, generated as random integers from 0 to 9 on the host and on the GPU. For each size it times (1) `gpuArray(host_data)` with `gputimeit`, (2) `gather(gpu_data)` with `gputimeit`, (3) `gpu_data+1.0` with `gputimeit`, and (4) `host_data+1.0` with `timeit`. Transfer bandwidth is size/time; the read/write estimates count two passes (2×size/time). The script converts these values to decimal GB/s and reports the maximum in each series. Its transfer plot distinguishes send from gather, using a logarithmic size axis and the default linear bandwidth axis. The read/write plot compares GPU and host arithmetic on logarithmic size and bandwidth axes.

The header records reference figures for an NVIDIA Tesla K40 (2013), Titan V (2017, TCC mode), and Tesla A100 (2021, PCIe) in a Dell PowerEdge T630: respectively send/gather/bw-GPU/bw-CPU values of 9.7/2.6/190/64, 10.3/2.6/568/64, and 10.4/2.6/1291/64. These are historical workstation notes, not guaranteed results; the comment does not give units for those four-tuples (the script's current console output labels its calculated bandwidths GB/s).

## Limits

Results depend on device, MATLAB/toolbox version, and host/device configuration. The function does not return the per-size timing or bandwidth arrays, and the source provides no device selection by name or input size override. Plotting uses Spinach helpers such as `kgrid`, `klegend`, `kfigure`, and `scale_figure`.
