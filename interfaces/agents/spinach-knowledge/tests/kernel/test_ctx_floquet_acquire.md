# tests/kernel/test_ctx_floquet_acquire.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_floquet_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_floquet_acquire.m)

## Purpose

Regression test for the Floquet context when used with `acquire()`. The test verifies that `floquet()` correctly projects states into Floquet space and runs `acquire()` on a tiny anisotropic one-spin MAS (magic-angle spinning) calculation, checking the returned time-domain trace for basic physical and dimensional invariants.

## Behaviour

1. Announces the test target with `fprintf('TESTING: Floquet acquire path\n')`.
2. Initialises a regression test result via `new_test_result('kernel/ctx_floquet_acquire', 'Floquet acquire path', 'floquet() must project states into Floquet space and run acquire().')`.
3. Builds a one-spin anisotropic Liouville-space spin system:
   - `sys.magnet = 14.1` (field strength in Tesla as given in the source).
   - `sys.isotopes = {'1H'}`.
   - `inter.zeeman.eigs = {[-2 -2 4]}` with `inter.zeeman.euler = {[0 0 0]}` (anisotropic shielding tensor).
   - `bas.formalism = 'sphten-liouv'`, `bas.approximation = 'none'`, `bas.projections = {+1}`.
   - The system is assembled with `test_spin_system(sys, inter, bas)`.
4. Sets up a tiny Floquet acquisition parameter struct:
   - `parameters.spins = {'1H'}`.
   - `parameters.rho0 = state(spin_system, 'L+', '1H')`.
   - `parameters.coil = state(spin_system, 'L+', '1H')`.
   - `parameters.decouple = {}`, `parameters.offset = 0`.
   - `parameters.sweep = 2000`, `parameters.npoints = 3`, `parameters.rate = 1000`.
   - `parameters.axis = [1 1 1]`, `parameters.max_rank = 1`.
   - `parameters.grid = 'leb_2ang_rank_5'`.
   - `parameters.serial = true`, `parameters.verbose = 0`.
5. Runs the production Floquet context: `fid = floquet(spin_system, @acquire, parameters, 'nmr')`.
6. Performs three checks:
   - **FID length:** `test_close(result, 'floquet FID length', numel(fid), parameters.npoints, 0, 0, ...)` — verifies `acquire()` returns one point per requested time sample.
   - **Zero-time signal:** computes `fid_zero = parameters.coil' * parameters.rho0` and checks `fid(1)` against it with tolerances `1e-12` (absolute and relative) — the Floquet-space projection must preserve the initial coil overlap.
   - **Finite FID:** `test_true(result, 'floquet finite FID', all(isfinite(fid(:))), ...)` — a short Floquet propagation should not produce NaN or Inf values.

## Inputs and outputs

**Syntax:**

```matlab
result = test_ctx_floquet_acquire()
```

**Inputs:** None. The function takes no arguments; all system and acquisition parameters are defined internally.

**Outputs:**

- `result` — regression test result object with explanatory messages, accumulated from the `test_close` and `test_true` checks described above.

## References

- Source file: [tests/kernel/test_ctx_floquet_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_floquet_acquire.m) in the Spinach repository.
