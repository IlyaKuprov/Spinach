# kernel/optimcon/hessreg.m

- Signature: `[H,data]=hessreg(spin_system,H,g,data)`
- Source: [kernel/optimcon/hessreg.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/hessreg.m)

## Purpose

Applies rational-function-optimisation (RFO) regularisation to a Newton-Raphson Hessian and gradient pair. It returns a regularised Hessian and the updated diagnostic counter; it does not return or modify the gradient. This helper does not perform a line search.

## Inputs and settings

- `H` must be a real, square, symmetric numeric matrix. `g` must be a real column vector with one element per Hessian dimension.
- The four RFO settings are read from `spin_system.control.reg_alpha`, `reg_phi`, `reg_max_iter`, and `reg_max_cond`. This function does not assign their defaults. The diagnostic structure must already provide `data.count.rfo`.

## Regularisation

If `H` is positive definite and its 2-norm condition number is already below `reg_max_cond`, the function returns it unchanged and takes no RFO iterations. Otherwise, for up to `reg_max_iter` iterations it forms the augmented matrix with blocks `alpha^2*H`, `alpha*g`, `alpha*g'`, and zero. It computes the smallest eigenvalue shift needed to make the augmented matrix nonnegative, subtracts that shift times the identity, removes the final row and column, and divides the remaining Hessian by `alpha^2`. It then multiplies `alpha` by `reg_phi`, increments `data.count.rfo`, and stops early if the condition number is below the target.

Finally, the result is replaced by its real symmetric part. If its condition number still meets or exceeds `reg_max_cond`, the routine emits a warning that the target was not reached. The source validates `H` and `g`, but does not validate or initialise the settings or counter structure.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=hessreg.m)
