# tests/kernel/test_operator_conversion_suite.m

## Purpose

Regression test for the Hilbert-to-Liouville operator conversion utilities in Spinach. It verifies direct vectorisation identities for left, right, commutation, and anticommutation superoperators, `unit_oper` dimensions in Zeeman Hilbert and Liouville formalisms, and `lindbladian` rate calibration.

## Behaviour

- Announces the test target with `fprintf('TESTING: Hilbert/Liouville conversion helpers\n')`.
- Initialises a `new_test_result` named `kernel/operator_conversion_suite` with the description "Hilbert/Liouville conversion helpers" and the requirement "operator conversion functions must implement vectorised product identities.".
- Builds a non-trivial Hermitian operator `H = S.z + 0.2*S.x` from `pauli(2)`, with `unit = speye(2)`.
- Checks `hilb2liouv` against direct Kronecker-product identities:
  - `'left'` against `kron(unit,H)`.
  - `'right'` against `kron(transpose(H),unit)`.
  - `'comm'` against `kron(unit,H) - kron(transpose(H),unit)`.
  - `'acomm'` against `kron(unit,H) + kron(transpose(H),unit)`.
  - `'statevec'` against `H(:)`.
  All comparisons use tolerances `1e-15` (absolute and relative) via `test_close`.
- Checks `unit_oper` dimensions using a test spin system with `sys.magnet=0`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`:
  - With `bas.formalism='zeeman-hilb'` and `bas.approximation='none'`, `unit_oper` is compared to `speye(2)`.
  - With `bas.formalism='zeeman-liouv'`, `unit_oper` is compared to `speye(4)`.
  Both use tolerances `1e-15`.
- Checks `lindbladian` calibration with `A_left=diag([1 0])`, `A_right=diag([0 1])`, `rho=[1;1]`, and `rate=3.5`. It computes `obs = real((rho'*R*rho)/(rho'*rho))` and compares it to `-rate` with tolerances `1e-12`, verifying that `lindbladian()` rescales the generator to the requested experimental decay rate.

## Inputs and outputs

- **Outputs**:
  - `result` — regression test result with explanatory messages, accumulated through repeated `test_close` calls.
- **Inputs**: None. The function takes no arguments.

## References

- Source: [tests/kernel/test_operator_conversion_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_operator_conversion_suite.m)
