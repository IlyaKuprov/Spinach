# tests/kernel/test_ctx_powder_average.m

## Purpose

Regression test that verifies the powder averaging path in `powder()`: the default powder average must equal an explicit weighted summation of the individual orientation traces using the spherical grid weights.

## Behaviour

- Announces the test target with `fprintf('TESTING: Powder weighted sum path\n')`.
- Initialises a regression test result via `new_test_result` with the identifier `kernel/ctx_powder_average`, the description `'Powder weighted sum path'`, and the specification `'powder() must sum orientation outputs with grid weights.'`.
- Builds a one-spin anisotropic Liouville-space spin system:
  - `sys.magnet = 14.1`
  - `sys.isotopes = {'1H'}`
  - `inter.zeeman.eigs = {[-2 -2 4]}`
  - `inter.zeeman.euler = {[0 0 0]}`
  - `bas.formalism = 'sphten-liouv'`
  - `bas.approximation = 'none'`
  - `bas.projections = {+1}`
  - The system is created with `test_spin_system(sys,inter,bas)`.
- Sets up a small powder acquisition with:
  - `parameters.spins = {'1H'}`
  - `parameters.rho0 = state(spin_system,'L+','1H')`
  - `parameters.coil = state(spin_system,'L+','1H')`
  - `parameters.decouple = {}`
  - `parameters.offset = 0`
  - `parameters.sweep = 2000`
  - `parameters.npoints = 3`
  - `parameters.grid = 'leb_2ang_rank_5'`
  - `parameters.serial = true`
  - `parameters.verbose = 0`
- Runs the averaged powder calculation: `[fid_avg,sph_grid] = powder(spin_system,@acquire,parameters,'nmr')`.
- Runs the per-orientation calculation by setting `parameters.sum_up = false`, returning `fid_cells` from `powder(spin_system,@acquire,parameters,'nmr')`.
- Assembles the reference weighted sum: `fid_ref = sph_grid.weights(1)*fid_cells{1}` accumulated over all orientations as `fid_ref = fid_ref + sph_grid.weights(n)*fid_cells{n}` for `n = 2:numel(fid_cells)`.
- Compares the averaged and reference results with `test_close(result,'weighted powder sum',fid_avg,fid_ref,1e-12,1e-12,...)`, with the message `'the powder average must equal the explicit grid-weighted sum of all orientations'`.

## Inputs and outputs

- **Inputs**: none. The function takes no arguments.
- **Outputs**: `result` — regression test result with explanatory messages, as stated in the file header.

## References

- Source file: [tests/kernel/test_ctx_powder_average.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_powder_average.m)
