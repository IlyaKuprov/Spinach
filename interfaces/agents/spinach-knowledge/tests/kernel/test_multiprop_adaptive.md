# tests/kernel/test_multiprop_adaptive.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_multiprop_adaptive.m`
- Signature: `result=test_multiprop_adaptive()`
- Total lines: 140

## Purpose

Tests adaptive repeated propagator application. Syntax: result=test_multiprop_adaptive()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Adaptive repeated propagator application\n')`.
- Lines 20-23: State the propagation target of the test; implemented by `result=new_test_result('kernel/multiprop_adaptive', 'Adaptive repeated propagator application', 'multiprop() must apply propagators repeatedly by binary adaptive squarin…`.
- Lines 25-26: Build the minimum Spinach system fields required by clean_up(); implemented by `spin_system.sys.enable={}`.
- Lines 33-34: Define a non-normal sparse propagator and a state vector; implemented by `P=sparse([0.90 0.10 0.00;0.00 0.80 0.20;0.05 0.00 0.70])`.
- Lines 38-39: Compare state-vector propagation with an explicit matrix power; implemented by `rho_obs=multiprop(spin_system,P,rho,N)`.
- Lines 44-45: Define a diagonal sparse propagator and a state vector; implemented by `P=spdiags([0.95;0.80;0.65],0,3,3)`.
- Lines 49-50: Compare sparse diagonal state-vector propagation with a matrix-power reference; implemented by `rho_obs=multiprop(spin_system,P,rho,N)`.
- Lines 55-56: Define a square Liouville-space vector stack; implemented by `P=sparse([0.90 0.10 0.00;0.00 0.80 0.20;0.05 0.00 0.70])`.
- Lines 60-61: Compare vector-stack propagation with an explicit matrix-power reference; implemented by `rho_obs=multiprop(spin_system,P,rho,N)`.
- Lines 66-67: Check the zero-step shortcut; implemented by `rho_obs=multiprop(spin_system,P,rho,0)`.
- Lines 71-72: Check scalar vector branch consistency; implemented by `spin_system.bas.formalism='zeeman-wavef'`.
- Lines 77-78: Define a Hilbert-space unitary propagator and a density matrix; implemented by `spin_system.bas.formalism='zeeman-hilb'`.
- Lines 85-86: Compare density-matrix propagation with an explicit matrix power; implemented by `P_ref=P^N`.
- Lines 92-93: Define a sparse non-unitary propagator and a sparse density matrix; implemented by `P=sparse([0.90 0.10 0.00;0.00 0.80 0.20;0.05 0.00 0.70])`.
- Lines 97-98: Compare sparse density-matrix propagation with an explicit matrix power; implemented by `P_ref=P^double(N)`.
- Lines 104-105: Define a diagonal sparse propagator and a density matrix; implemented by `P=spdiags(exp(1i*[0.10;0.20;0.30]),0,3,3)`.
- Lines 109-110: Compare sparse diagonal density-matrix propagation with an explicit matrix power; implemented by `P_ref=P^double(N)`.
- Lines 116-117: Enable clean-up after propagator squaring; implemented by `spin_system.bas.formalism='zeeman-liouv'`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/multiprop_adaptive', 'Adaptive repeated propagator application', 'multiprop() must apply propagators repeatedly by binary adaptive squarin…`.
- Lines 26: computes `spin_system.sys.enable` using `spin_system.sys.enable={}`.
- Lines 27: computes `spin_system.sys.disable` using `spin_system.sys.disable={}`.
- Lines 28: computes `spin_system.bas.formalism` using `spin_system.bas.formalism='zeeman-liouv'`.
- Lines 29: computes `spin_system.tols.prop_chop` using `spin_system.tols.prop_chop=0`.
- Lines 30: computes `spin_system.tols.dense_matrix` using `spin_system.tols.dense_matrix=0.15`.
- Lines 31: computes `spin_system.tols.small_matrix` using `spin_system.tols.small_matrix=200`.
- Lines 34: computes `P` using `P=sparse([0.90 0.10 0.00;0.00 0.80 0.20;0.05 0.00 0.70])`.
- Lines 35: computes `rho` using `rho=[1.0;2.0;3.0]`.
- Lines 36: computes `N` using `N=uint64(13)`.
- Lines 39: computes `rho_obs` using `rho_obs=multiprop(spin_system,P,rho,N)`.
- Lines 40: computes `rho_ref` using `rho_ref=(P^double(N))*rho`.
- Lines 79: computes `S` using `S=pauli(2)`.
- Lines 80: computes `H` using `H=2*pi*(0.30*S.x+0.70*S.z)`.
- Lines 86: computes `P_ref` using `P_ref=P^N`.
- Lines 129: computes `caught` using `caught=false`.

## Outputs

- result -regression test result with explanatory messages
- The test checks binary adaptive squaring in multiprop() against explicit
- matrix-power references for state vectors and density matrices, and
- verifies clean-up after propagator squaring.

## Implementation structure

- Tests adaptive repeated propagator application. Syntax:
- result=test_multiprop_adaptive()
- result -regression test result with explanatory messages
- The test checks binary adaptive squaring in multiprop() against explicit
- matrix-power references for state vectors and density matrices, and
- verifies clean-up after propagator squaring.
- Announce the test target
- State the propagation target of the test
- Build the minimum Spinach system fields required by clean_up()
- Define a non-normal sparse propagator and a state vector
- Compare state-vector propagation with an explicit matrix power
- Define a diagonal sparse propagator and a state vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `multiprop()`, `uint64()`, `double()`, `test_close()`, `spdiags()`, `pauli()`, `speye()`, `test_true()`.
