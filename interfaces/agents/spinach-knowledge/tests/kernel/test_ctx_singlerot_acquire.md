# tests/kernel/test_ctx_singlerot_acquire.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_singlerot_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_singlerot_acquire.m)

## Purpose

Regression test for the single-rotor context (`singlerot()`) when used with `acquire()`. It runs a tiny anisotropic one-spin MAS calculation and checks the returned time-domain trace for basic physical and dimensional invariants.

## Behaviour

- Announces the test target with `fprintf('TESTING: Single-rotor acquire path\n')`.
- Initialises a regression test result via `new_test_result('kernel/ctx_singlerot_acquire', ...)`, stating that `singlerot()` must project states into rotor space and run `acquire()`.
- Builds a one-spin anisotropic Liouville-space system:
  - `sys.magnet=14.1`, `sys.isotopes={'1H'}`.
  - `inter.zeeman.eigs={[-2 -2 4]}`, `inter.zeeman.euler={[0 0 0]}`.
  - `bas.formalism='sphten-liouv'`, `bas.approximation='none'`, `bas.projections={+1}`.
  - The system is constructed with `test_spin_system(sys,inter,bas)`.
- Sets up a tiny MAS acquisition with:
  - `parameters.spins={'1H'}`.
  - `parameters.rho0=state(spin_system,'L+','1H')` and `parameters.coil=state(spin_system,'L+','1H')`.
  - `parameters.decouple={}`, `parameters.offset=0`.
  - `parameters.sweep=2000`, `parameters.npoints=3`.
  - `parameters.rate=1000`, `parameters.axis=[1 1 1]`.
  - `parameters.max_rank=1`, `parameters.grid='single_crystal'`.
  - `parameters.serial=true`, `parameters.verbose=0`.
- Runs the production single-rotor context: `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Performs three checks:
  1. FID length: `numel(fid)` must equal `parameters.npoints` (tolerances 0 and 0), because `acquire()` should return one point for each requested time sample.
  2. Zero-time signal: `fid(1)` must match `parameters.coil'*parameters.rho0` (tolerances `1e-12` and `1e-12`), because rotor-space projection must preserve the initial coil overlap.
  3. Finite FID: `all(isfinite(fid(:)))` must be true, because short MAS propagation should not produce NaN or Inf values.

## Inputs and outputs

```matlab
result=test_ctx_singlerot_acquire()
```

- **Outputs:** `result` — regression test result with explanatory messages.
- **Inputs:** none.

## References

- Source file: [tests/kernel/test_ctx_singlerot_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_singlerot_acquire.m)
