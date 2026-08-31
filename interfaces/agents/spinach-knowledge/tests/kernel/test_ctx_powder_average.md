# tests/kernel/test_ctx_powder_average.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_powder_average.m`
- Signature: `result=test_ctx_powder_average()`
- Total lines: 65

## Purpose

Tests powder averaging against explicit weighted summation. Syntax: result=test_ctx_powder_average()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Powder weighted sum path\n')`.
- Lines 19-22: State the powder-averaging target of the test; implemented by `result=new_test_result('kernel/ctx_powder_average', 'Powder weighted sum path', 'powder() must sum orientation outputs with grid weights.')`.
- Lines 24-25: Build a one-spin anisotropic Liouville-space system; implemented by `sys.magnet=14.1`.
- Lines 34-35: Set up a tiny powder acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 46-47: Run the averaged powder calculation; implemented by `[fid_avg,sph_grid]=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 49-50: Run the per-orientation powder calculation; implemented by `parameters.sum_up=false`.
- Lines 53-54: Assemble the independent weighted sum; implemented by `fid_ref=sph_grid.weights(1)*fid_cells{1}`.
- Lines 59-61: Check the powder average accumulation path; implemented by `result=test_close(result,'weighted powder sum',fid_avg,fid_ref,1e-12,1e-12, 'the powder average must equal the explicit grid-weighted sum of all orientations')`.

### Control flow inferred from the code

- Line 55: `for` loop over `n=2:numel(fid_cells)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/ctx_powder_average', 'Powder weighted sum path', 'powder() must sum orientation outputs with grid weights.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-2 -2 4]}`.
- Lines 28: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]}`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `bas.projections` using `bas.projections=+1`.
- Lines 32: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 36: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 37: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 38: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 39: computes `parameters.offset` using `parameters.offset=0`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=2000`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=3`.
- Lines 42: computes `parameters.grid` using `parameters.grid='leb_2ang_rank_5'`.
- Lines 43: computes `parameters.serial` using `parameters.serial=true`.

## Outputs

- result -regression test result with explanatory messages
- The test asks powder() for individual orientation traces and checks that
- the default powder average is the same weighted sum.

## Implementation structure

- Tests powder averaging against explicit weighted summation. Syntax:
- result=test_ctx_powder_average()
- result -regression test result with explanatory messages
- The test asks powder() for individual orientation traces and checks that
- the default powder average is the same weighted sum.
- Announce the test target
- State the powder-averaging target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a tiny powder acquisition
- Run the averaged powder calculation
- Run the per-orientation powder calculation
- Assemble the independent weighted sum

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `powder()`, `test_spin_system()`, `state()`, `test_close()`.
