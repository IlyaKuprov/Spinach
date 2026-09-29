# tests/kernel/test_liquid_single_spin_fid.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_liquid_single_spin_fid.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_liquid_single_spin_fid.m)

## Purpose

Regression test for a one-spin liquid-state free induction decay (FID). The test verifies that, with no Hamiltonian and no relaxation, transverse magnetisation is constant in time — i.e. a zero-offset isolated spin has a constant FID.

## Behaviour

- Announces the test target with `fprintf('TESTING: Single-spin liquid-state FID\n')`.
- Registers a new test result via `new_test_result` under the name `kernel/liquid_single_spin_fid`, described as "Single-spin liquid-state FID", with the property that "a zero-offset isolated spin has a constant free induction decay."
- Builds a one-spin Liouville-space spin system:
  - `sys.magnet=14.1` (14.1 T magnet), `sys.isotopes={'1H'}`.
  - `inter.zeeman.scalar={0}` (zero Zeeman scalar coupling).
  - `bas.formalism='sphten-liouv'`, `bas.approximation='none'`.
  - The system is constructed with `test_spin_system(sys,inter,bas)`.
- Sets up a zero-offset acquisition:
  - `parameters.spins={'1H'}`.
  - `parameters.rho0=state(spin_system,'L+','1H')` and `parameters.coil=state(spin_system,'L+','1H')`.
  - `parameters.decouple={}`, `parameters.offset=0`.
  - `parameters.sweep=1000`, `parameters.npoints=8`, `parameters.zerofill=8`.
  - `parameters.axis_units='Hz'`, `parameters.invert_axis=0`.
- Simulates the FID with `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Closes the test with `test_close(result,'constant zero-offset FID',fid,fid(1)*ones(size(fid)),1e-12,1e-12,...)`, checking that the simulated signal equals its first point repeated across all points, to absolute and relative tolerances of `1e-12`, with the explanation "without precession or relaxation the detected coherence is time-independent".

## Inputs and outputs

**Syntax:**

```matlab
result=test_liquid_single_spin_fid()
```

**Outputs:**

- `result` — regression test result structure with explanatory messages, as produced by `new_test_result` and finalised by `test_close`.

**Inputs:**

- None.

## References

- Spinach kernel test: `tests/kernel/test_liquid_single_spin_fid.m` ([GitHub source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_liquid_single_spin_fid.m)).
- Related functions used by this test: `new_test_result`, `test_spin_system`, `state`, `liquid`, `acquire`, `test_close`.
