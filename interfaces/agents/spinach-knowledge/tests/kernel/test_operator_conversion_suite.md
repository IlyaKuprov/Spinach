# tests/kernel/test_operator_conversion_suite.m

- Signature: `result=test_operator_conversion_suite()`

## Purpose

Tests Hilbert-to-Liouville operator conversion identities, unit-operator dimensions, and Lindbladian rate calibration.

## Physical / mathematical content

For the Hermitian two-level operator `H=S.z+0.2*S.x`, checks that left and right multiplication vectorise as `kron(I,H)` and `kron(transpose(H),I)`. Commutation and anticommutation use their difference and sum, respectively; state-vector conversion stacks columns as `H(:)`.

## Numerical / algorithmic content

Compares all five `hilb2liouv` modes with their direct expressions using absolute and relative tolerances of `1e-15`. Checks that `unit_oper` is `speye(2)` in `zeeman-hilb` and `speye(4)` in `zeeman-liouv` for a single `1H` spin. With `A_left=diag([1 0])`, `A_right=diag([0 1])`, `rho=[1;1]`, and `rate=3.5`, checks that the real expectation value of `lindbladian(A_left,A_right,rho,rate)` is `-rate` within `1e-12`.

## Outputs

- `result` — regression test result with explanatory messages.