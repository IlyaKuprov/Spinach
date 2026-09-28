# kernel/optimcon/hessreg.m

- Signature: `[H,data]=hessreg(spin_system,H,g,data)`

## Purpose

Regularises a real symmetric Newton-Raphson Hessian using rational function optimisation (RFO), shifting its spectrum as needed to obtain a better-conditioned Hessian.

## Parameters / inputs

- `spin_system` — Spinach system object; regularisation settings are read from `spin_system.control`.
- `H` — real symmetric Hessian matrix.
- `g` — real column gradient with the same number of elements as the dimension of `H`.
- `data` — diagnostic structure whose `data.count.rfo` field is incremented for each RFO iteration.

## Outputs

- `H` — regularised Hessian. If it is already positive definite and below the configured condition-number limit, the input Hessian is returned unchanged.
- `data` — diagnostic structure with the RFO iteration count updated.

## Implementation

The routine uses `reg_alpha`, `reg_phi`, `reg_max_iter`, and `reg_max_cond` from `spin_system.control`. Each iteration forms the augmented RFO Hessian, shifts by its lowest eigenvalue when needed, and checks the resulting Hessian's condition number. It symmetrises the final result and warns if the target condition number was not reached.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=hessreg.m)
