# experiments/pseudocon/hfc2pms.m

- Signature: `[pms,pms_tensor]=hfc2pms(A,chi,isotope)`
- MATLAB source: [`experiments/pseudocon/hfc2pms.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/hfc2pms.m)

## Purpose

Calculates the full paramagnetic-shift tensor, including contact and pseudocontact contributions, from the hyperfine-coupling and magnetic-susceptibility tensors. The relation is Equation 10 in [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).

## Inputs and outputs

- `A` — real symmetric 3-by-3 hyperfine-coupling tensor in Gauss, normalised per unpaired electron in the `S*A*I` spin-Hamiltonian convention (the convention returned by `gparse.m`).
- `chi` — real symmetric 3-by-3 magnetic-susceptibility tensor in cubic ångströms (Å³).
- `isotope` — isotope label as a character string, for example `'1H'`; its nuclear gyromagnetic ratio is obtained with `spin(isotope)`.
- `pms_tensor` — paramagnetic-shift tensor in ppm.
- `pms` — isotropic part, calculated as one third of the tensor trace, in ppm.

The source notes that Gauss is used for the hyperfine coupling because this convention does not depend on the electron g-tensor.

## Calculation

The code uses the full supplied tensors in `pms_tensor = 1e6*(1/(4*pi))*A*chi/C`, with `pms = trace(pms_tensor)/3`. Unlike the PCS-only calculation, this routine does not first discard isotropic components, so its reported paramagnetic shift includes the contact as well as pseudocontact contribution.

The code obtains `gamma_n=spin(isotope)`, and uses `hbar=1.05457173e-34` and `mu0=4*pi*1e-7`. Its conversion factor is `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.

The input checks require both tensors to be real symmetric 3-by-3 matrices and `isotope` to be a character string. The routine transforms tensors; it is not a pulse-sequence function.

## References

- [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G)
- [Spin Dynamics Wiki: hfc2pms.m](https://spindynamics.org/wiki/index.php?title=hfc2pms.m)
