# kernel/utilities/xyz2hfc.m

## Purpose

Converts point electron and nuclear coordinates into a hyperfine interaction tensor.

## Behaviour

- Syntax: `A=xyz2hfc(exyz,nxyz,isotope)`.
- The function first runs a consistency check (`grumble`) on the inputs.
- Fundamental constants used: `hbar=1.054571730e-34` and `mu0=4*pi*1e-7`.
- The nuclear magnetogyric ratio is obtained via `spin(isotope)`.
- The electron position is used as the origin: `nxyz=nxyz-exyz`.
- A prefactor is computed as `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.
- The dipolar matrix is `D=3*(nxyz'*nxyz)/norm(nxyz,2)^5-eye(3)/norm(nxyz,2)^3`.
- The returned tensor is `A=C*D`.
- Gauss units are used for hyperfine couplings because they do not depend on the electron g-tensor.
- The tensor returned is the one that enters the spin Hamiltonian as `S*A*I`; it does not scale with the number of unpaired electrons because the electron spin operator already carries that magnitude.
- Input validation errors:
  - `exyz must be a 1x3 real row vector.`
  - `nxyz must be a 1x3 real row vector.`
  - `e_xyz and n_xyz coordinates must be different.` (raised when `norm(nxyz-exyz,2)==0`)
  - `isotope specification must be a character string.`

## Inputs and outputs

Inputs:

- `exyz` — Cartesian coordinates of the electron, a 1x3 row vector in Angstrom.
- `nxyz` — Cartesian coordinates of the nucleus, a 1x3 row vector in Angstrom.
- `isotope` — isotope specification, e.g. `'13C'`.

Output:

- `A` — hyperfine coupling tensor, Gauss.

## References

- Source: [kernel/utilities/xyz2hfc.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/xyz2hfc.m)
- Wiki: [xyz2hfc.m](https://spindynamics.org/wiki/index.php?title=xyz2hfc.m)
