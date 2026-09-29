# experiments/pseudocon/xyz2pms.m

- Source: [experiments/pseudocon/xyz2pms.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/xyz2pms.m)
- Wiki: [xyz2pms.m](https://spindynamics.org/wiki/index.php?title=xyz2pms.m)
- Signature: `pms_tensor=xyz2pms(nxyz,sxyz,chi)`

## Purpose

Computes the paramagnetic shielding tensor generated at one nuclear coordinate by a point magnetic-susceptibility centre. This is a geometric forward calculation, not a shielding-tensor fit.

## Inputs and coordinate convention

- `nxyz` and `sxyz` are real numeric row vectors `[x y z]` for the nucleus and susceptibility centre, respectively, in Angstroms. The displacement row is `r=nxyz-sxyz`.
- `chi` is the susceptibility tensor in cubic Angstroms, supplied as a real `3×3` matrix or five components ordered `[chi(1,1), chi(1,2), chi(1,3), chi(2,2), chi(2,3)]`. The five-component form is expanded to a symmetric traceless matrix, using `chi(3,3)=-chi(1,1)-chi(2,2)`.

## Calculation and output

The code forms `D=3*(r'*r)/norm(r,2)^5-eye(3)/norm(r,2)^3` and returns `pms_tensor=1e6*(1/(4*pi))*D*chi`. The output is a `3×3` paramagnetic shielding tensor in ppm. The matrix order is specifically `D*chi`; the source notes this order follows Spinach's convention of placing the magnetic field on the left in the Zeeman Hamiltonian. With a full matrix input, that matrix is used as supplied; the validator does not enforce symmetry or tracelessness. There is no coincident-site guard, and the formula is singular at `r=0`.
