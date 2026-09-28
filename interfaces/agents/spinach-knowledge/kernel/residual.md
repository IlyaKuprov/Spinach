# kernel/residual.m

- Signature: `spin_system=residual(spin_system)`

## Purpose

Replaces anisotropic Zeeman and spin-spin coupling tensors with their isotropic parts plus weak residual-order contributions, using the user-supplied order matrix for each chemical subsystem.

## Physical / mathematical content

For a tensor `T`, the function separates its isotropic part `iso` and evaluates the residual contribution as `extra_zz=trace(order_matrix*(T-iso))`. It then forms `iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3])` for the updated tensor.

## Numerical / algorithmic content

The function checks that `spin_system.inter.order_matrix` is present and nonempty, then processes Zeeman and coupling tensors within each chemical subsystem. A coupling tensor whose matrix 2-norm is below `2*pi*spin_system.tols.inter_cutoff` after replacement is removed.

## Parameters / inputs

- `spin_system` - output of `create.m` containing the spin system and interaction data, including the order matrix.

## Outputs

- `spin_system` - the input structure with its interaction tensors overwritten by their partial-order residual forms.

## Notes

- Applicable to weak residual order in high-field NMR spectroscopy.
- Any required relaxation superoperator must be computed before this function, because the interaction tensors are overwritten.
- The `liquid.m` context invokes this function automatically when `parameters.needs` contains `'rdc'`.

## Source

- ledwards@cbs.mpg.de
- ilya.kuprov@weizmann.ac.il
- <https://spindynamics.org/wiki/index.php?title=residual.m>
