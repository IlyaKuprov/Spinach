# tests/kernel/test_ctx_gridfree_acquire.m

## Purpose

Regression test for the grid-free Fokker-Planck context, verifying that `gridfree()` correctly projects spin system states into SLE (stochastic Liouville equation) space and runs `acquire()` on the resulting context.

## Behaviour

The test announces its target with `fprintf('TESTING: Grid-free acquire path\n')` and registers a test result via `new_test_result()` with the identifier `kernel/ctx_gridfree_acquire`, description `'Grid-free acquire path'`, and the requirement that `gridfree()` must project states into SLE space and run `acquire()`.

It then builds a one-spin anisotropic Liouville-space system using `test_spin_system()` with:

- `sys.magnet = 14.1`
- `sys.isotopes = {'1H'}`
- `inter.zeeman.eigs = {[-2 -2 4]}` (anisotropic Zeeman eigenvalues)
- `inter.zeeman.euler = {[0 0 0]}`
- `bas.formalism = 'sphten-liouv'`
- `bas.approximation = 'none'`
- `bas.projections = {+1}`

A tiny grid-free acquisition is configured with `parameters` containing `spins = {'1H'}`, `rho0` and `coil` both set to `state(spin_system,'L+','1H')`, `decouple = {}`, `offset = 0`, `sweep = 2000`, `npoints = 3`, `rate = 1000`, `axis = [1 1 1]`, `max_rank = 2`, and `verbose = 0`.

The production grid-free context is invoked as:

```matlab
fid = gridfree(spin_system,@acquire,parameters,'nmr');
```

Three checks are then performed:

1. **FID length**: `numel(fid)` is compared against `parameters.npoints` with zero tolerances, requiring `acquire()` to return one point per requested time sample.
2. **Zero-time signal**: `fid(1)` is compared against `fid_zero = parameters.coil'*parameters.rho0` with absolute and relative tolerances of `1e-12`, verifying that SLE-space projection preserves the initial coil overlap.
3. **Finite FID**: `all(isfinite(fid(:)))` is checked, ensuring short grid-free propagation produces no NaN or Inf values.

## Inputs and outputs

```matlab
result = test_ctx_gridfree_acquire()
```

- **Output**: `result` — regression test result object with explanatory messages, accumulated through `test_close()` and `test_true()` checks.
- **Input**: none.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_gridfree_acquire.m)
