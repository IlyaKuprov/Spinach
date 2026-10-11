# Spinach regression tests

This directory contains physically motivated regression tests for the Spinach kernel and its core numerical machinery. Each test file states the physical or mathematical property under test, how the Spinach calculation is constructed, why the reference answer is correct, and what numerical tolerance is acceptable. Reference answers are analytic, algebraic, or derived explicitly in the test file.

Run from the Spinach root:

```matlab
addpath('tests');
results=run_tests();
```

Useful variants:

```matlab
list_tests();
run_test('kernel/pauli_spin_half_algebra');
run_tests('pattern','relaxation');
run_tests('verbose',true);
```


The CUDA sparse-product comparison is registered in the main manifest:

```matlab
run_test('kernel/cuda_product');
```

It compares the gateway directly with independent CPU products for random
square and rectangular sparse matrices in all real/complex combinations.
GPU checks report explicit `SKIP` messages when a usable GPU or the compiled
gateway is unavailable, following the existing GPU coverage convention.
With both prerequisites present, numerical and non-loader CUDA errors fail
the test. Build the gateway with `compile_mex()` before running GPU coverage.

GPU-dependent custom CUDA sparse-product regressions (not part of the CPU manifest):

```matlab
addpath('kernel/overloads','tests/kernel');
test_alg3_fallback();
test_cuda_sparse();
```

These require a supported GPU. The fallback test uses an isolated wrapper copy
against missing and invalid binaries and mocked non-loader failures; shipped
binaries are never modified. The product test requires the compiled gateway
and checks independent CPU references, real/complex combinations, empty and
rectangular inputs, cancellation, bounded-memory row splitting, input
immutability, and ordinary MATLAB operations on the returned sparse gpuArray.
