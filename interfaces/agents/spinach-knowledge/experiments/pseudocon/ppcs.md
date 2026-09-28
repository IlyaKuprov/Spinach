# experiments/pseudocon/ppcs.m

- Signature: `pcs=ppcs(nxyz,sxyz,chi)`

## Purpose

Computes the pseudocontact shift, in ppm, at each nuclear coordinate in `nxyz` for a point susceptibility centre at `sxyz`.

## Physical / mathematical content

The susceptibility tensor `chi` is supplied as a real `3x3` matrix or as five components ordered `[chi(1,1) chi(1,2) chi(1,3) chi(2,2) chi(2,3)]`. For five components, the routine constructs a symmetric traceless tensor, with `chi(3,3)=-chi(1,1)-chi(2,2)`. Coordinates are taken relative to `sxyz`; the PCS is evaluated with an `l=2` spherical-harmonic expansion proportional to `r^-3` and converted to ppm.

## Parameters / inputs

- `nxyz` — real nuclear coordinates, one `[x y z]` row per nucleus, in Angstroms.
- `sxyz` — real `[x y z]` susceptibility-centre coordinate in Angstroms.
- `chi` — real `3x3` susceptibility tensor in cubic Angstroms, or the five components listed above.

## Output

- `pcs` — predicted pseudocontact shift in ppm at each nucleus.
