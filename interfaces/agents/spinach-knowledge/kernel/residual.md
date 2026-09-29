# kernel/residual.m

- Signature: `spin_system=residual(spin_system)`

## Purpose

Replaces the Zeeman and spin-spin coupling interaction tensors with partial-order residual forms using the supplied order matrix. The source describes this for weak residual order in high-field NMR spectroscopy in a liquid crystal.

## Physical / mathematical content

For each chemical subsystem, the corresponding order-matrix entry is used to compute a residual scalar from the traceless part of a tensor: `extra_zz=trace(order_matrix{s}*(T-iso))`, where `iso=trace(T)*eye(3)/3`. The replacement is `iso+diag([-extra_zz/3, -extra_zz/3, 2*extra_zz/3])`. This operation is applied to the subsystem's Zeeman tensors and pair-coupling tensors.

## Numerical / algorithmic content

The order matrix must be present and nonempty. Coupling tensors whose matrix 2-norm, after replacement, is below `2*pi*spin_system.tols.inter_cutoff` are cleared. The function overwrites the interaction tensors in `spin_system`; compute any required relaxation superoperator before calling it. The `liquid.m` context calls it automatically when `parameters.needs` contains `'rdc'`.

## Parameters / inputs

- `spin_system` — output of `create.m` with interaction data and the order matrix. Adjustable parameters are set during `create.m`.

## Outputs

- `spin_system` — the same structure with the affected interaction tensors replaced by their partial-order residual forms.

## Source

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/residual.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=residual.m)
