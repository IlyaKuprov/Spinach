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


GPU-dependent SpGEMM gateway compatibility regression (not part of the CPU manifest):

```matlab
addpath('kernel/overloads','tests/kernel');
test_alg3_fallback();
```

This requires a supported GPU and tests an isolated production-wrapper copy
against missing and invalid MEX binaries and mocked non-loader failures.
Shipped binaries are never renamed, modified, or rebuilt.
