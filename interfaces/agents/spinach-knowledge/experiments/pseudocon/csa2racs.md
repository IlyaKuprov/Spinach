# experiments/pseudocon/csa2racs.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/csa2racs.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=csa2racs.m)

## Purpose

Computes a high-temperature estimate of the residual anisotropic chemical shift from supplied chemical-shift-anisotropy and magnetic-susceptibility tensors. The source identifies the expression with Equation (2) of [Otting et al.](https://doi.org/10.1021/ja0564259). This is a tensor calculation, not a pulse-transfer or MAS simulation.

## Inputs

- `csa` is a real numeric 3-by-3 chemical-shift tensor in ppm.
- `chi` is a real numeric 3-by-3 magnetic-susceptibility tensor in cubic Angstroms.
- `B` is a real scalar magnetic induction in Tesla.
- `T` is a positive real scalar absolute temperature in Kelvin.

The routine keeps the spherical-rank-2 component of each tensor, contracts them by a trace, and applies the field- and temperature-dependent factor. In the source this is `1e-30*(B^2/(15*mu_0*k_b*T))*trace(csa*chi)`, with `mu_0=4*pi*1e-7` and `k_b=1.38064852e-23`. The scalar field enters as `B^2`; no field-orientation or pulse-delay input appears in this signature.

## Output

- `racs` is the residual anisotropic chemical shift in ppm.
