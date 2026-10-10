# kernel/overloads/cuda_sparse_by_sparse.m

Source: [kernel/overloads/cuda_sparse_by_sparse.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/cuda_sparse_by_sparse.m)

## Purpose and interface

`C=cuda_sparse_by_sparse(A,B,alg)` computes a sparse GPU matrix product using cuSPARSE SpGEMM. Both operands must be sparse double gpuArrays with compatible dimensions. Real and complex operands may be mixed; the result is a sparse double gpuArray, complex if either input is complex and neither is all-zero, as in native `mtimes`.

`alg` is a required CPU double scalar: 1, 2, or 3 selects `CUSPARSE_SPGEMM_ALG1`, `ALG2`, or `ALG3`. `propagator` uses ALG2, which was the fastest on R2026b and H200 for menthol propagator squarings (5822-dimensional, 1.7M and 4.4M nonzeros); ALG3 took about twice as long.

## Numerical representation

R2026b stores a sparse gpuArray as zero-based CSR of the matrix itself with int32 indices. `cuda_sparse_by_sparse_mex` reads the value, column, and row-offset device arrays in place from this undocumented layout, never writes them, and validates it before use: all three must be device allocations on the current GPU spanning the matrix, the first row offset must be zero, and the last must equal `nzmax`, the stored entry count (GPU arithmetic can leave explicit zeros, so `nnz` may be smaller).

ALG1 and ALG2 run on MATLAB's 32-bit indices. ALG3 gets 64-bit index copies because 32-bit ALG3 corrupts device memory on large products; its chunk fraction is halved from one until its estimation and compute buffers fit into free device memory. When cuSPARSE or CUDA runs out of resources, the rows of A are bisected at half of their nonzeros; a row block is passed as pointer offsets into MATLAB's arrays with a rebased copy of its row offsets. A real operand of a mixed product is copied into complex values. The product returns as one-based int32 row-major triplets, assembled with GPU `sparse`, which drops explicit zeros. Products with more than `intmax('int32')` nonzeros are rejected.

This explicit utility does not override native `gpuArray` multiplication. Dense, CPU, or incompatible operands are rejected rather than silently converted.

## Compiled dependency

`cuda_sparse_by_sparse_mex` must be compiled and on the MATLAB path to use cuSPARSE. The shipped Linux `.mexa64` targets MATLAB R2026b with CUDA 13 and has been tested on H200; the binary depends on `libcusparse.so.12`, `libcudart.so.13`, and `libcuda.so.1`. The wrapper uses native GPU multiplication when the platform MEX is absent, when its invocation raises `MATLAB:mex:ErrInvalidMEXFile`, or when the gateway raises `Spinach:cuda_sparse_by_sparse_mex:layout` because the storage layout of the running MATLAB release is not recognised. All input validation remains active, and other failures, including CUDA computation and allocation errors, propagate unchanged.
