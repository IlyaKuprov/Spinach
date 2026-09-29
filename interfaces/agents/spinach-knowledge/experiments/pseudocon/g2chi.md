# experiments/pseudocon/g2chi.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/g2chi.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=g2chi.m)

## Purpose

Calculates a high-temperature Curie-law magnetic-susceptibility tensor from an input g tensor, absolute temperature, and electron spin quantum number. It does not fit an NMR transfer experiment or use a magnetic-field input.

## Inputs

- `g` is a real numeric 3-by-3 dimensionless g tensor; the Curie-law prefactor already contains the square of the Bohr magneton, so do not multiply `g` by that constant.
- `T` is a positive real scalar temperature in Kelvin.
- `S` is a positive real scalar integer or half-integer spin. The source examples are 1/2, 1, and 3/2.

The code forms `prefactor = S*(S+1)*mu_0*mu_b^2/(3*k_b*T)` and returns `chi = 1e30*prefactor*(g*transpose(g))`. Its constants are `mu_b=9.274009994e-24`, `mu_0=4*pi*1e-7`, and `k_b=1.38064852e-23`.

## Output

- `chi` is a 3-by-3 magnetic-susceptibility tensor in cubic Angstroms.
