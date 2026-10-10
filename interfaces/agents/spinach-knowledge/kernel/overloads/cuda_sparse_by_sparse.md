# kernel/overloads/cuda_sparse_by_sparse.m

Source: [kernel/overloads/cuda_sparse_by_sparse.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/cuda_sparse_by_sparse.m)

## Purpose and interface

`C=cuda_sparse_by_sparse(A,B)` multiplies compatible sparse double gpuArrays using custom CUDA CSR arithmetic, without cuSPARSE SpGEMM. All four real/complex combinations are supported without promoting or copying real input values. The result is a MATLAB-owned sparse double gpuArray, ready for ordinary MATLAB operations; neither input is modified. The obsolete cuSPARSE algorithm-selector argument has been removed.

## Algorithm and storage

Gustavson row-wise multiplication separates structural counting from numerical accumulation. Narrow column ranges use a bitset for the pattern and a dense shared-memory accumulator for values. Each warp owns a disjoint column interval in the numeric phase, so no atomics are needed there. Wider ranges use shared-memory hash tables sized from the row's structural upper bound or measured output count. Overflowing symbolic tasks bisect their column range; no intermediate scalar-product array or dense matrix is allocated. Shared-memory limits derive from the device properties, and global task metadata scale with rows and output structure.

The gateway reads R2026b's undocumented zero-based CSR buffers in place. It validates device ownership, allocation extents, row boundaries, and strictly increasing column indices. GPU arithmetic may retain stored zeros and unused allocation capacity. The active structural entry count comes from the final CSR row offset; `nzmax` bounds the allocation capacity, while `nnz` counts numerical nonzeros. These quantities need not coincide. Storage uses MATLAB's int32 ABI; loop indices, task ranges, output sums, and pointer arithmetic use int64. Outputs exceeding the int32 storage limit are rejected rather than truncated.

A fresh sparse gpuArray is allocated from the symbolic GPU pattern with placeholder values, using MATLAB's `sparse` constructor. The numeric kernel writes directly into that object's owned CSR values; no numerical result is copied or reformatted. The constructor still needs temporary pattern buffers. Structural entries whose numerical sums cancel to zero may remain stored; MATLAB `nnz`, `gather`, and Spinach `clean_up` handle them.

## Availability and errors

The shipped Linux binary targets R2026b and H200 (sm_90), with CUDA runtime and driver dependencies but no direct cuSPARSE dependency. Missing or unloadable platform binaries retain native GPU multiplication. An unrecognised internal layout raises `Spinach:cuda_sparse_by_sparse_mex:layout`; this does not fall back to native multiplication. CUDA computation and allocation errors propagate.

Algorithm references: Gustavson (1978), DOI [10.1145/355791.355796](https://doi.org/10.1145/355791.355796); Davis et al., [Sparse Direct Methods, Algorithm 2.2](https://people.engr.tamu.edu/davis/publications_files/survey_tech_report.pdf); [GraphBLAS saxpy3](https://github.com/DrTimothyAldenDavis/GraphBLAS/blob/stable/Source/mxm/GB_AxB_saxpy3.c); [spECK](https://github.com/GPUPeople/spECK).
