# tests/kernel/test_ctx_doublerot_acquire.m

## Purpose

Regression test for the double-rotor context in Spinach, exercising the `acquire()` simulation method through `doublerot()`. The test runs a tiny anisotropic one-spin double-rotation calculation and checks the returned time-domain trace for basic physical and dimensional invariants.

Source: [tests/kernel/test_ctx_doublerot_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_doublerot_acquire.m)

## Behavior

- Announces the test target with `fprintf('TESTING: Double-rotor acquire path\n')`.
- Initializes a regression test result via `new_test_result('kernel/ctx_doublerot_acquire', ...)`, stating that `doublerot()` must project states into double-rotor space and run `acquire()`.
- Builds a one-spin anisotropic Liouville-space spin system with:
  - `sys.magnet = 14.1`
  - `sys.isotopes = {'1H'}`
  - Zeeman interaction eigenvalues `[-2 -2 4]` with Euler angles `[0 0 0]`
  - `bas.formalism = 'sphten-liouv'`, `bas.approximation = 'none'`, `bas.projections = {+1}`
- Sets up a double-rotation acquisition with:
  - `parameters.rho0` and `parameters.coil` both set to `state(spin_system,'L+','1H')`
  - `parameters.offset = 0`, `parameters.sweep = 2000`, `parameters.npoints = 3`
  - `parameters.rate_outer = 800`, `parameters.rate_inner = 2400`
  - `parameters.rank_outer = 1`, `parameters.rank_inner = 1`
  - `parameters.axis_outer = [sqrt(2/3) 0 sqrt(1/3)]`
  - `parameters.axis_inner = [sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30)]`
  - `parameters.grid = 'single_crystal'`, `parameters.serial = true`, `parameters.verbose = 0`
- Runs the production double-rotor context: `fid = doublerot(spin_system,@acquire,parameters,'nmr')`.
- Performs three checks:
  1. **FID length**: `numel(fid)` equals `parameters.npoints` (tolerances 0 and 0), because `acquire()` should return one point per requested time sample.
  2. **Zero-time signal**: `fid(1)` matches `fid_zero = parameters.coil'*parameters.rho0` (tolerances `1e-12` and `1e-12`), because double-rotor projection must preserve the initial coil overlap.
  3. **Finite FID**: `all(isfinite(fid(:)))` is true, because short double-rotor propagation should not produce NaN or Inf values.

## Inputs and outputs

```matlab
result = test_ctx_doublerot_acquire()
```

- **Output**: `result` — regression test result structure with explanatory messages, accumulated through `test_close` and `test_true` checks.
- **Input**: none.

## References

- Source file: [tests/kernel/test_ctx_doublerot_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_doublerot_acquire.m) (Spinach repository, main branch).
