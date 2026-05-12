# Spinach regression tests

This directory contains hand-written, physically motivated regression tests for the Spinach kernel and its core numerical machinery.  The tests are intentionally readable: each test file states the physical or mathematical property under test, how the Spinach calculation is constructed, why the reference answer is correct, and what numerical tolerance is acceptable.

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

The suite avoids opaque generated wrappers and opaque `.mat` known-answer blobs.  Reference answers are analytic, algebraic, or derived explicitly in the test file.
