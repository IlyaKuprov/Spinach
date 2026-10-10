# tests/kernel/test_cuda_sparse.m

Source: [tests/kernel/test_cuda_sparse.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cuda_sparse.m)

`test_cuda_sparse()` exercises the production custom CUDA sparse multiplication gateway on a supported GPU. Independent CPU products provide exact integer references and normwise floating-point references for all four real/complex operand combinations. The test includes rectangular and empty shapes, cancellation and explicitly stored zeros, a million-column sparse hash, and a dense row that requires column-range splitting. Inputs must remain unchanged.

The returned sparse gpuArray is also checked through transpose, absolute value, matrix-vector multiplication, triplet extraction, repeated custom multiplication, and Spinach's GPU cleanup. Run it separately as documented in `tests/README.md`; it is not part of the CPU-only manifest.
