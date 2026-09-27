# experiments/pseudocon/hfc2pcs.m

- Signature: `[pcs,pcs_tensor]=hfc2pcs(A,chi,isotope)`

## Purpose

Converts a hyperfine coupling tensor and magnetic susceptibility tensor to the pseudocontact shift tensor and its isotropic part, excluding the contact contribution, according to Equation 10 of [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).

## Parameters / inputs

- `A` — real symmetric 3-by-3 hyperfine coupling tensor in Gauss, normalised per unpaired electron in the `S*A*I` convention (as returned by `gparse.m`).
- `chi` — real symmetric 3-by-3 magnetic susceptibility tensor in Å^3.
- `isotope` — isotope label as a character string, for example `'1H'`.

## Outputs

- `pcs_tensor` — pseudocontact shift tensor in ppm.
- `pcs` — isotropic pseudocontact shift in ppm, calculated as `trace(pcs_tensor)/3`.

## Method

The routine keeps only the rank-2 components of `A` and `chi`, obtains the nuclear gyromagnetic ratio from `spin(isotope)`, and evaluates the full shift tensor using the fundamental-constant factor in Equation 10. It then takes one third of the tensor trace for the isotropic shift. Inputs are checked for real symmetric 3-by-3 tensors and a character-string isotope label.

## References

- [10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G)
- [Spin Dynamics Wiki: hfc2pcs.m](https://spindynamics.org/wiki/index.php?title=hfc2pcs.m)
