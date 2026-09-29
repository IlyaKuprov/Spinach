# tests/kernel/test_ctx_doublerot_acquire.m

## Purpose

Regression test for the double-rotor context in Spinach, exercising the `acquire()` simulation method through `doublerot()`. The test runs a tiny anisotropic one-spin double-rotation calculation and checks the returned time-domain trace for basic physical and dimensional invariants.

Source: [tests/kernel/test_ctx_doublerot_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_doublerot_acquire.m)

## Test configuration and invariants

A single `1H` spin at 14.1 T with axial Zeeman eigenvalues `[-2 -2 4]` is represented without basis truncation in `sphten-liouv`. The single-crystal context uses outer/inner rotor rates of 800/2400 Hz, rank one for each rotor, and axes `[sqrt(2/3) 0 sqrt(1/3)]` and `[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))]`. With a 2000 Hz sweep it acquires three time points.

The test checks three properties of double-rotor `acquire()` propagation: the FID has exactly the requested three points; its first point equals the initial coil–state overlap within `1e-12` absolute and relative tolerance; and all returned points are finite. It does not compare a full spectrum with experimental data.

## Inputs and outputs

```matlab
result = test_ctx_doublerot_acquire()
```

- **Output**: `result` — regression test result structure with explanatory messages, accumulated through `test_close` and `test_true` checks.
- **Input**: none.

## References

- Source file: [tests/kernel/test_ctx_doublerot_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_doublerot_acquire.m) (Spinach repository, main branch).
