# experiments/pseudocon/ppcs.m

- Source: [experiments/pseudocon/ppcs.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/ppcs.m)
- Wiki: [ppcs.m](https://spindynamics.org/wiki/index.php?title=ppcs.m)
- Signature: `pcs=ppcs(nxyz,sxyz,chi)`

## Purpose

Calculates the point-centre pseudocontact shift at one or more nuclear coordinates. This is a forward dipolar-shift calculation, not a tensor-fitting routine.

## Inputs and coordinate convention

- `nxyz` is a real numeric `N×3` array of nuclear coordinates `[x y z]`, in Angstroms.
- `sxyz` is one real numeric centre coordinate `[x y z]`, in Angstroms. The displacement used is `nxyz-sxyz`.
- `chi` is the magnetic-susceptibility tensor in cubic Angstroms: either a real `3×3` matrix or five real components ordered `[chi(1,1), chi(1,2), chi(1,3), chi(2,2), chi(2,3)]`. The five-value form is expanded to a symmetric traceless matrix with `chi(3,3)=-chi(1,1)-chi(2,2)`.

## Calculation and output

The routine converts each displacement to spherical coordinates, obtains rank-2 irreducible components with `qform2sph`, and accumulates the spherical-harmonic expansion `(1/(4*pi))*chi2(-m+3)*r^(-3)*spher_harmon(2,m,theta,phi)` for `m=[2,1,0,-1,-2]`. It returns `1e6*real(sum)`, a column of predicted PCS values in ppm, one per row of `nxyz`. Thus the spatial dependence is inverse-cubic and the calculated shift comes from the rank-2 susceptibility contribution.

The input checks require real numeric coordinates with three columns for `nxyz`, a real numeric three-element row for `sxyz`, and a real `3×3` tensor after any five-component expansion. No coincident-site check is made; the inverse-cubic expression is singular when a nuclear coordinate equals the centre.
