# tests/kernel/test_ctx_powder_crystal.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_powder_crystal.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ctx_powder_crystal.m)

## Purpose

Regression test that verifies the equivalence of the `powder()` and `crystal()` static simulation contexts at a single orientation. The test uses the `single_crystal` grid so that `powder()` and `crystal()` represent the same Euler orientation of the same anisotropic one-spin Hamiltonian.

## Behaviour

- Announces the test target with `fprintf('TESTING: Powder and crystal single orientation\n')`.
- Initialises a regression test result via `new_test_result('kernel/ctx_powder_crystal', ...)`, stating that `powder()` with the `single_crystal` grid must match `crystal()`.
- Builds a one-spin anisotropic Liouville-space spin system with `test_spin_system(sys,inter,bas)` using:
  - `sys.magnet=14.1`, `sys.isotopes={'1H'}`.
  - `inter.zeeman.eigs={[-2 -2 4]}`, `inter.zeeman.euler={[0 0 0]}`.
  - `bas.formalism='sphten-liouv'`, `bas.approximation='none'`, `bas.projections={+1}`.
- Sets up a short static acquisition with `parameters.spins={'1H'}`, `rho0` and `coil` both set to `state(spin_system,'L+','1H')`, `decouple={}`, `offset=0`, `sweep=2000`, `npoints=4`, `grid='single_crystal'`, `orientation=[0 0 0]`, `serial=true`, `verbose=0`.
- Runs both production static contexts: `fid_powder=powder(spin_system,@acquire,parameters,'nmr')` and `fid_crystal=crystal(spin_system,@acquire,parameters,'nmr')`.
- Compares the two trajectories with `test_close(result,'single orientation powder',fid_powder,fid_crystal,1e-12,1e-12,...)`, with the explanation that the `single_crystal` powder grid has unit weight at zero Euler angles.

## Inputs and outputs

```matlab
result=test_ctx_powder_crystal()
```

- **Output**: `result` — regression test result with explanatory messages.
- **Input**: none.

## References
