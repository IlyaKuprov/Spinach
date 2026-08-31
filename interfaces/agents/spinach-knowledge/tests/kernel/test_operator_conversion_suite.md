# tests/kernel/test_operator_conversion_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_operator_conversion_suite.m`
- Signature: `result=test_operator_conversion_suite()`
- Total lines: 67

## Purpose

Tests Hilbert-to-Liouville operator conversion utilities. Syntax: result=test_operator_conversion_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Hilbert/Liouville conversion helpers\n')`.
- Lines 20-23: State the operator-conversion target of the test; implemented by `result=new_test_result('kernel/operator_conversion_suite', 'Hilbert/Liouville conversion helpers', 'operator conversion functions must implement vectorised product ident…`.
- Lines 25-26: Define a non-trivial Hermitian operator; implemented by `S=pauli(2)`.
- Lines 30-32: Check direct vectorisation formulas; implemented by `result=test_close(result,'hilb2liouv left',hilb2liouv(H,'left'),kron(unit,H),1e-15,1e-15, 'left multiplication vectorises as kron(I,H)')`.
- Lines 42-43: Check unit operator dimensions in major formalisms; implemented by `sys.magnet=0`.
- Lines 56-57: Check Lindbladian calibration on a simple vector; implemented by `A_left=diag([1 0])`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/operator_conversion_suite', 'Hilbert/Liouville conversion helpers', 'operator conversion functions must implement vectorised product ident…`.
- Lines 26: computes `S` using `S=pauli(2)`.
- Lines 27: computes `H` using `H=S.z+0.2*S.x`.
- Lines 28: computes `unit` using `unit=speye(2)`.
- Lines 43: computes `sys.magnet` using `sys.magnet=0`.
- Lines 44: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 45: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 46: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 47: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 48: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 57: computes `A_left` using `A_left=diag([1 0])`.
- Lines 58: computes `A_right` using `A_right=diag([0 1])`.
- Lines 59: computes `rho` using `rho=[1;1]`.
- Lines 60: computes `rate` using `rate=3.5`.
- Lines 61: computes `R` using `R=lindbladian(A_left,A_right,rho,rate)`.
- Lines 62: computes `obs` using `obs=real((rho'*R*rho)/(rho'*rho))`.

## Outputs

- result -regression test result with explanatory messages
- The test checks direct vectorisation identities for left, right,
- commutation, and anticommutation superoperators, unit_oper dimensions,
- and Lindbladian rate calibration.

## Implementation structure

- Tests Hilbert-to-Liouville operator conversion utilities. Syntax:
- result=test_operator_conversion_suite()
- result -regression test result with explanatory messages
- The test checks direct vectorisation identities for left, right,
- commutation, and anticommutation superoperators, unit_oper dimensions,
- and Lindbladian rate calibration.
- Announce the test target
- State the operator-conversion target of the test
- Define a non-trivial Hermitian operator
- Check direct vectorisation formulas
- Check unit operator dimensions in major formalisms
- Check Lindbladian calibration on a simple vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `pauli()`, `speye()`, `test_close()`, `hilb2liouv()`, `transpose()`, `test_spin_system()`, `unit_oper()`, `lindbladian()`.
