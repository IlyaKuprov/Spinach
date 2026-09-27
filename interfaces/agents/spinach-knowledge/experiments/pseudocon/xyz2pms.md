# experiments/pseudocon/xyz2pms.m

- Signature: `pms_tensor=xyz2pms(nxyz,sxyz,chi)`

## Purpose

Computes the paramagnetic shielding tensor, in ppm, at a nuclear coordinate `nxyz` due to a point magnetic-susceptibility centre at `sxyz`.

## Physical / mathematical content

The displacement from the centre is `r=nxyz-sxyz`. The routine forms the dipolar matrix `D=3*r'*r/|r|^5-I/|r|^3` and returns `pms_tensor=10^6*D*chi/(4*pi)`. The full susceptibility tensor is used, including its isotropic part; the matrix order `D*chi` follows the Spinach Zeeman-Hamiltonian convention noted in the source.

## Parameters / inputs

- `nxyz` — real nuclear coordinate `[x y z]` in Angstroms (a row vector).
- `sxyz` — real susceptibility-centre coordinate `[x y z]` in Angstroms.
- `chi` — magnetic susceptibility tensor in cubic Angstroms, supplied as a `3x3` matrix or five components ordered `[chi(1,1) chi(1,2) chi(1,3) chi(2,2) chi(2,3)]`. With five components, the routine constructs a symmetric traceless tensor.

## Output

- `pms_tensor` — predicted paramagnetic shielding tensor in ppm.
