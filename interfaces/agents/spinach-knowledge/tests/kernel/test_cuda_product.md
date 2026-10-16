# CUDA sparse-product regression

`test_cuda_product()` returns a standard Spinach regression result and is
registered as `kernel/cuda_product` in the main unit test suite. It compares
`cuda_sparse_by_sparse_mex` directly with CPU multiplication of the gathered
operands, rather than the wrapper that can fall back to native multiplication.
Three square/rectangular shapes cover all four real/complex operand pairings,
with independently generated imaginary components and a reproducible random
stream that is restored on exit.

The absolute tolerance is a complex dot-product roundoff bound scaled by
`norm(abs(A)*abs(B),'fro')`, accounting for both CPU and GPU accumulation.
Checks also require native sparse GPU output. A MATLAB release before R2026a, missing usable GPU, compute capability below 7.0, or missing gateway
produces an explicit skip message; the specific
`MATLAB:mex:ErrInvalidMEXFile` loader error is also skipped. Other CUDA errors
and incorrect products remain failures.

Run `run_test('kernel/cuda_product')` after `compile_mex()` to exercise it.
