# tests/kernel/test_dynamic_state_projection_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_state_projection_suite.m`
- Signature: `result=test_dynamic_state_projection_suite()`
- Total lines: 110

## Purpose

Tests dynamic state projection helper paths. Syntax: result=test_dynamic_state_projection_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Dynamic state projection helpers\n')`.
- Lines 19-22: State the state-projection target of the test; implemented by `result=new_test_result('kernel/dynamic_state_projection_suite', 'Dynamic state projection helpers', 'state projection helpers must preserve Hermitian adjoint pairs, stat…`.
- Lines 24-25: Build a two-deuteron Hilbert-space system; implemented by `sys.magnet=0`.
- Lines 33-34: Request all deuteron-pair coherences; implemented by `[~,~,~,Tc,Qc]=deut_pair(spin_d,1,2)`.
- Lines 36-38: Check adjoint pairing of triplet coherences; implemented by `result=test_close(result,'triplet coherence adjoint 1',Tc{1},Tc{3}',1e-12,1e-12, 'T0->T- and T-->T0 coherences must be Hermitian adjoints')`.
- Lines 42-44: Check adjoint pairing of quintet coherences; implemented by `result=test_close(result,'quintet coherence adjoint 1',Qc{1},Qc{5}',1e-12,1e-12, 'Q-->Q-- and Q--->Q- coherences must be Hermitian adjoints')`.
- Lines 52-53: Request dephased deuteron-pair populations; implemented by `options.dephasing=1`.
- Lines 56-57: Build the longitudinal two-spin Hamiltonian used by the dephasing claim; implemented by `Hz=operator(spin_d,'Lz',1)+operator(spin_d,'Lz',2)`.
- Lines 59-61: Check that dephased population states are stationary under Az+Bz; implemented by `result=test_close(result,'dephased singlet stationarity',Hz*Sd-Sd*Hz,0*Sd,1e-12,1e-12, 'dephased deuteron-pair singlet population must commute with Az+Bz')`.
- Lines 71-72: Build a spherical-tensor system for captured stateinfo() output; implemented by `sys_s.magnet=0`.
- Lines 81-82: Capture and inspect the printed state composition report; implemented by `printed=evalc('stateinfo(spin_s,state(spin_s,''Lz'',''1H''),1);')`.
- Lines 90-91: Build a triplet-electron Hilbert-space system; implemented by `sys_e.magnet=0`.
- Lines 99-100: Project an isotropic zero-field population through arbitrary high-field axes; implemented by `ZFS=[2 0.1 0.2; 0.1 1 0.3; 0.2 0.3 -3]*1e6`.
- Lines 104-106: Check that an isotropic triplet population remains maximally mixed; implemented by `result=test_close(result,'isotropic zftrip projection',rho_zf,speye(3)/3,1e-12,1e-12, 'an isotropic zero-field triplet population is invariant under projection and depha…`.

### Control flow inferred from the code

- Line 62: `for` loop over `n=1:numel(Td)`.
- Line 66: `for` loop over `n=1:numel(Qd)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_state_projection_suite', 'Dynamic state projection helpers', 'state projection helpers must preserve Hermitian adjoint pairs, stat…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'2H','2H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0,0}`.
- Lines 28: computes `inter.temperature` using `inter.temperature=300`.
- Lines 29: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_d` using `spin_d=test_spin_system(sys,inter,bas)`.
- Lines 34: computes `[~,~,~,Tc,Qc]` using `[~,~,~,Tc,Qc]=deut_pair(spin_d,1,2)`.
- Lines 53: computes `options.dephasing` using `options.dephasing=1`.
- Lines 54: computes `[Sd,Td,Qd]` using `[Sd,Td,Qd]=deut_pair(spin_d,1,2,options)`.
- Lines 57: computes `Hz` using `Hz=operator(spin_d,'Lz',1)+operator(spin_d,'Lz',2)`.
- Lines 72: computes `sys_s.magnet` using `sys_s.magnet=0`.
- Lines 73: computes `sys_s.isotopes` using `sys_s.isotopes={'1H'}`.
- Lines 74: computes `inter_s.zeeman.scalar` using `inter_s.zeeman.scalar={0}`.
- Lines 75: computes `inter_s.temperature` using `inter_s.temperature=300`.
- Lines 76: computes `bas_s.formalism` using `bas_s.formalism='sphten-liouv'`.
- Lines 77: computes `bas_s.approximation` using `bas_s.approximation='none'`.

## Outputs

- result -regression test result with explanatory messages
- The test checks deuteron-pair coherences, dephased population stationarity,
- captured stateinfo() output, and isotropic zero-field triplet projection.

## Implementation structure

- Tests dynamic state projection helper paths. Syntax:
- result=test_dynamic_state_projection_suite()
- result -regression test result with explanatory messages
- The test checks deuteron-pair coherences, dephased population stationarity,
- captured stateinfo() output, and isotropic zero-field triplet projection.
- Announce the test target
- State the state-projection target of the test
- Build a two-deuteron Hilbert-space system
- Request all deuteron-pair coherences
- Check adjoint pairing of triplet coherences
- Check adjoint pairing of quintet coherences
- Request dephased deuteron-pair populations

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `deut_pair()`, `test_close()`, `operator()`, `int2str()`, `evalc()`, `stateinfo()`, `state()`, `test_true()`, `contains()`, `report()`, `zftrip()`, `speye()`.
