# experiments/pseudocon/hfc2pcs.m

- Signature: `[pcs,pcs_tensor]=hfc2pcs(A,chi,isotope)`
- MATLAB source: [`experiments/pseudocon/hfc2pcs.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/hfc2pcs.m)

## Purpose

Computes the pseudocontact-shift (PCS) tensor from hyperfine-coupling and magnetic-susceptibility tensors, excluding the isotropic contact contribution. The relation is the PCS part of Equation 10 in [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).

## Inputs and outputs

- `A` — real symmetric 3-by-3 hyperfine-coupling tensor in Gauss, normalised per unpaired electron in the `S*A*I` spin-Hamiltonian convention (the convention returned by `gparse.m`).
- `chi` — real symmetric 3-by-3 magnetic-susceptibility tensor in cubic ångströms (Å³).
- `isotope` — isotope label as a character string, for example `'1H'`; its nuclear gyromagnetic ratio is obtained with `spin(isotope)`.
- `pcs_tensor` — 3-by-3 PCS tensor in ppm.
- `pcs` — isotropic part, calculated as one third of the tensor trace, in ppm.

The source notes that Gauss is used for the hyperfine coupling because this convention does not depend on the electron g-tensor.

## Calculation

Before multiplying the tensors, the function converts each to spherical-tensor components and reconstructs it from rank 2 only. Thus the calculation uses the anisotropic components of both inputs to form the PCS tensor and excludes the contact contribution. With `A2` and `chi2` denoting those reconstructed rank-2 tensors, the implemented calculation is `pcs_tensor = 1e6*(1/(4*pi))*A2*chi2/C`, followed by `pcs = trace(pcs_tensor)/3`.

The code obtains `gamma_n=spin(isotope)`, and uses `hbar=1.05457173e-34` and `mu0=4*pi*1e-7`. Its conversion factor is `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.

The input checks require both tensors to be real symmetric 3-by-3 matrices and `isotope` to be a character string. The routine performs a tensor conversion; it is not a pulse-sequence function.

## References

- [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G)
- [Spin Dynamics Wiki: hfc2pcs.m](https://spindynamics.org/wiki/index.php?title=hfc2pcs.m)
