# tests/kernel/test_dynamic_iserstep_highorder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_iserstep_highorder.m`
- Signature: `result=test_dynamic_iserstep_highorder()`
- Total lines: 89

## Purpose

Tests nonlinear high-order iserstep branches. Syntax: result=test_dynamic_iserstep_highorder()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `local_hilb_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Nonlinear high-order iserstep branches\n')`.
- Lines 20-23: State the nonlinear Lie-step target of the test; implemented by `result=new_test_result('kernel/dynamic_iserstep_highorder', 'Nonlinear high-order iserstep branches', 'high-order Lie and RKMK solvers must execute nonlinear generators…`.
- Lines 25-26: Build a one-proton Hilbert-space spin system; implemented by `spin_system=local_hilb_system()`.
- Lines 32-33: Build a Hermitian density matrix with non-zero coherences; implemented by `rho=[0.7,0.2+0.1i;0.2-0.1i,0.3]`.
- Lines 35-36: Define a mildly nonlinear, non-commuting Hamiltonian field; implemented by `H0=0.31*Sx+0.17*Sz`.
- Lines 42-43: Check the explicit zero-time shortcut in LG4A; implemented by `rho_zero=iserstep(spin_system,{Lfun,0,'LG4A'},rho,0)`.
- Lines 47-48: Build a refined DP8 reference using two half steps; implemented by `rho_ref=iserstep(spin_system,{Lfun,0,'RKMK-DP8'},rho,dt/2)`.
- Lines 51-52: Exercise nonlinear high-order branches against the refined reference; implemented by `methods={'LG4A','RKMK4','RKMK-DP5','RKMK-DP8','RKMK-RKF45'}`.
- Lines 56-57: Take one nonlinear Lie step with this method; implemented by `rho_obs=iserstep(spin_system,{Lfun,0,methods{n}},rho,dt)`.
- Lines 59-61: Compare to the refined high-order reference; implemented by `result=test_close(result,['iserstep nonlinear ' methods{n}],rho_obs,rho_ref,tolerances(n),tolerances(n), 'nonlinear high-order branches agree with a refined DP8 referenc…`.
- Lines 63-65: Check density-matrix invariants preserved by unitary propagation; implemented by `result=test_close(result,['iserstep trace ' methods{n}],trace(rho_obs),trace(rho),1e-12,1e-12, 'unitary Hamiltonian propagation preserves the density-matrix trace')`.

### Control flow inferred from the code

- Line 54: `for` loop over `n=1:numel(methods)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_iserstep_highorder', 'Nonlinear high-order iserstep branches', 'high-order Lie and RKMK solvers must execute nonlinear generators…`.
- Lines 26: computes `spin_system` using `spin_system=local_hilb_system()`.
- Lines 27: computes `S` using `S=pauli(2)`.
- Lines 28: computes `Sx` using `Sx=S.x`.
- Lines 29: computes `Sy` using `Sy=S.y`.
- Lines 30: computes `Sz` using `Sz=S.z`.
- Lines 33: computes `rho` using `rho=[0.7,0.2+0.1i;0.2-0.1i,0.3]`.
- Lines 36: computes `H0` using `H0=0.31*Sx+0.17*Sz`.
- Lines 37: computes `H1` using `H1=-0.23*Sy+0.05*Sx`.
- Lines 38: computes `H2` using `H2=0.11*Sx-0.07*Sz`.
- Lines 39: computes `Lfun` using `Lfun=@(t,state)H0+sin(4*t)*H1+0.2*real(trace(Sz*state))*H2`.
- Lines 40: computes `dt` using `dt=5e-3`.
- Lines 43: computes `rho_zero` using `rho_zero=iserstep(spin_system,{Lfun,0,'LG4A'},rho,0)`.
- Lines 48: computes `rho_ref` using `rho_ref=iserstep(spin_system,{Lfun,0,'RKMK-DP8'},rho,dt/2)`.
- Lines 52: computes `methods` using `methods={'LG4A','RKMK4','RKMK-DP5','RKMK-DP8','RKMK-RKF45'}`.
- Lines 53: computes `tolerances` using `tolerances=[5e-8 5e-8 5e-9 5e-10 5e-9]`.
- Lines 57: computes `rho_obs` using `rho_obs=iserstep(spin_system,{Lfun,0,methods{n}},rho,dt)`.

### Local helper functions

- Line 73: `local_hilb_system()` — `function spin_system=local_hilb_system()`. Specify the spin system
  - Representative operation: `sys.magnet=0`.
  - Representative operation: `sys.isotopes={'1H'}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks zero-step handling, nonlinear generator execution, and
- agreement of high-order Lie and RKMK branches against a refined DP8
- reference on a compact Hilbert-space problem.

## Implementation structure

- Tests nonlinear high-order iserstep branches. Syntax:
- result=test_dynamic_iserstep_highorder()
- result -regression test result with explanatory messages
- The test checks zero-step handling, nonlinear generator execution, and
- agreement of high-order Lie and RKMK branches against a refined DP8
- reference on a compact Hilbert-space problem.
- Announce the test target
- State the nonlinear Lie-step target of the test
- Build a one-proton Hilbert-space spin system
- Build a Hermitian density matrix with non-zero coherences
- Define a mildly nonlinear, non-commuting Hamiltonian field
- Check the explicit zero-time shortcut in LG4A

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_hilb_system()`, `pauli()`, `iserstep()`, `test_close()`, `tolerances()`, `test_spin_system()`.
