# kernel/overloads/cuda_sparse_by_sparse.m

Source: [kernel/overloads/cuda_sparse_by_sparse.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/cuda_sparse_by_sparse.m)

## Purpose and interface

`C=cuda_sparse_by_sparse(A,B,chunk_fraction)` computes a sparse GPU matrix product using cuSPARSE `CUSPARSE_SPGEMM_ALG3`. Both operands must be sparse double gpuArrays with compatible dimensions. Real and complex operands may be mixed; the result is a sparse double gpuArray, complex whenever either input is complex. Matrix dimensions and input nonzero counts must fit into int32, and the gateway rejects an output nonzero count exceeding int32.

`chunk_fraction` is a required finite real CPU double scalar in `[realmin('single'),1]`. It controls ALG3 intermediate-product chunking, not a strict bound on total device memory. `propagator` uses `0.02` at its sparse GPU Taylor and squaring call sites.

## Numerical representation

MATLAB `find` exposes GPU-resident triplets in column order. The CUDA helper interprets this ordering as CSR storage of the transposed operands, computes `B.'*A.'`, and returns one-based triplets for reconstruction with GPU `sparse`. The transpose is not a conjugate transpose. Numerical array data stay on the GPU; dimensions and the chunk fraction are CPU metadata. Empty products retain the requested shape and complexity without running SpGEMM.

This explicit utility does not override native `gpuArray` multiplication. Dense, CPU, or incompatible operands are rejected rather than silently converted.

## Compiled dependency

`cuda_sparse_by_sparse_mex` must be compiled and on the MATLAB path, including for empty products. The shipped Linux `.mexa64` targets MATLAB R2026b with CUDA 13.1 and has been tested on H200. Compatible builds are required for other platforms or releases; the binary depends on `libcusparse.so.12` and `libcudart.so.13`. The CUDA source uses the public MATLAB GPU API and cuSPARSE. The wrapper reports a missing compiled helper; it does not silently fall back to native multiplication.
