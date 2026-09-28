# experiments/pseudocon/hfc2pms.m

- Signature: `[pms,pms_tensor]=hfc2pms(A,chi,isotope)`

## Purpose

Calculates the paramagnetic shift tensor, including contact and pseudocontact contributions, from the hyperfine coupling and magnetic susceptibility tensors, following Equation 10 of [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).

## Parameters / inputs

- `A` — real symmetric 3-by-3 hyperfine coupling tensor in Gauss, normalised per unpaired electron in the `S*A*I` convention (as returned by `gparse.m`). Gauss units avoid dependence on the electron g-tensor.
- `chi` — real symmetric 3-by-3 magnetic susceptibility tensor in Å^3.
- `isotope` — isotope label as a character string, for example `'1H'`.

## Outputs

- `pms_tensor` — paramagnetic shift tensor in ppm.
- `pms` — isotropic paramagnetic shift in ppm, calculated as `trace(pms_tensor)/3`.

## Method

The routine obtains the nuclear gyromagnetic ratio from `spin(isotope)` and evaluates the full tensor as `1e6*(1/(4*pi))*A*chi/C`, with `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`, using the constants defined in the source. It computes the isotropic shift from one third of the tensor trace. Inputs are checked for real symmetric 3-by-3 tensors and a character-string isotope label.

## References

- [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G)
- [Spin Dynamics Wiki: hfc2pms.m](https://spindynamics.org/wiki/index.php?title=hfc2pms.m)
