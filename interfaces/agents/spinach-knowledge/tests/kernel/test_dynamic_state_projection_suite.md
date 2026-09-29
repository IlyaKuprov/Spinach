# tests/kernel/test_dynamic_state_projection_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_state_projection_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_state_projection_suite.m)

## Purpose

Regression test suite for the dynamic state projection helpers in Spinach. The suite verifies that state projection helpers preserve Hermitian adjoint pairs, stationarity, and normalisation, covering deuteron-pair coherences, dephased population stationarity, captured `stateinfo()` output, and isotropic zero-field triplet projection.

## Behaviour

- Announces the test target with `fprintf('TESTING: Dynamic state projection helpers\n')` and initialises a regression test result via `new_test_result` for `kernel/dynamic_state_projection_suite`.
- Builds a two-deuteron Hilbert-space system with `sys.magnet=0`, `sys.isotopes={'2H','2H'}`, zero Zeeman scalars, `inter.temperature=300`, `bas.formalism='zeeman-hilb'`, and `bas.approximation='none'`, using `test_spin_system`.
- Requests all deuteron-pair coherences with `[~,~,~,Tc,Qc]=deut_pair(spin_d,1,2)` and checks adjoint pairing of triplet coherences (`Tc{1}` vs `Tc{3}'`, `Tc{2}` vs `Tc{4}'`) and quintet coherences (`Qc{1}` vs `Qc{5}'`, `Qc{2}` vs `Qc{6}'`, `Qc{3}` vs `Qc{7}'`, `Qc{4}` vs `Qc{8}'`), each with tolerances `1e-12`.
- Requests dephased deuteron-pair populations with `options.dephasing=1` and `[Sd,Td,Qd]=deut_pair(spin_d,1,2,options)`; builds the longitudinal two-spin Hamiltonian `Hz=operator(spin_d,'Lz',1)+operator(spin_d,'Lz',2)` and checks that the dephased singlet, each triplet, and each quintet population commutes with `Az+Bz` (tolerances `1e-12`).
- Builds a spherical-tensor single-proton system (`sys_s.isotopes={'1H'}`, `bas_s.formalism='sphten-liouv'`, `bas_s.approximation='none'`, `inter_s.temperature=300`, zero Zeeman scalar, `sys_s.magnet=0`) with `spin_s.sys.output=1`; captures the printed state composition report via `evalc('stateinfo(spin_s,state(spin_s,''Lz'',''1H''),1);')` and asserts the captured text contains `'state vector 2-norm'`, `'most populated basis states'`, and `'+5.000e-01'`.
- Builds a triplet-electron Hilbert-space system (`sys_e.isotopes={'E3'}`, `bas_e.formalism='zeeman-hilb'`, `bas_e.approximation='none'`, `inter_e.temperature=300`, zero Zeeman scalar, `sys_e.magnet=0`); projects an isotropic zero-field population through arbitrary high-field axes using `ZFS=[2 0.1 0.2; 0.1 1 0.3; 0.2 0.3 -3]*1e6`, `Z=diag([27.8 28.1 28.4])*1e9`, and `rho_zf=zftrip(spin_e,ZFS,[1/3 1/3 1/3],Z,0.05,1)`; checks that the isotropic triplet population remains maximally mixed by comparing `rho_zf` to `speye(3)/3` with tolerances `1e-12`.

## Inputs and outputs

Syntax:

```matlab
result=test_dynamic_state_projection_suite()
```

- **Output:** `result` — regression test result with explanatory messages, accumulated through `test_close` and `test_true` checks.
- **Input:** none.

## References

- [Spinach GitHub repository — test_dynamic_state_projection_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_state_projection_suite.m)
