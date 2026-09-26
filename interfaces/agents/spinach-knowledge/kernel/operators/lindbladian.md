# kernel/operators/lindbladian.m

- Signature: `R=lindbladian(A_left,A_right,rho,rlx_rate)`

## Purpose

Generates a Lindblad superoperator from user-specified left-side and right-side product superoperators and calibrates it using the experi- mental relaxation rate of a user-specified state. Syntax: R=lindbladian(A_left,A_right,rho,rlx_rate)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- A_left -left side product superoperator of the
- interaction that is causing relaxation
- (see operator.m and hamiltonian.m)
- A_right -right side product superoperator of the
- same interaction
- rho -the state vector whose relaxation rate
- is known from the experiment
- rlx_rate -experimental relaxation rate of rho

## Outputs

- R -Lindblad relaxation superoperator indu-
- ced by the interaction A, such that
- <rho|R|rho>/norm(rho,2)^2 = -rlx_rate

## Implementation structure

- Generates a Lindblad superoperator from user-specified left-side and
- right-side product superoperators and calibrates it using the experi-
- mental relaxation rate of a user-specified state. Syntax:
- R=lindbladian(A_left,A_right,rho,rlx_rate)
- A_left -left side product superoperator of the
- interaction that is causing relaxation
- (see operator.m and hamiltonian.m)
- A_right -right side product superoperator of the
- same interaction
- rho -the state vector whose relaxation rate
- is known from the experiment
- rlx_rate -experimental relaxation rate of rho
