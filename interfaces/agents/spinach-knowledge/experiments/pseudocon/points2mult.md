# experiments/pseudocon/points2mult.m

- Signature: `Ilm=points2mult(xyz,mxyz,rho,L,method)`

## Purpose

Computes multipole moments from a set of points with user-specified spin populations. Syntax: Ilm=points2mult(xyz,mxyz,rho,L,method) The multipoles in question are described in Equation 32 of:

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

## Parameters / inputs

- xyz -coordinates as [x y z] with multiple rows,
- at which density rho is evaluated, in Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z], in
- Angstroms.
- L -array of ranks of spherical harmonics of the
- probability density
- rho -column of the densities at the points xyz
- method -'points' or 'grid'. If the spin density is supplied
- as Mulliken spin populations at individual nuclei,
- choose 'points'; if the spin density is supplied as
- a probability on a uniform cubic grid obtained from
- ndgrid() function and vectorised, use 'grid'.
- Output:
- Ilm -multipole moments of the probability density

## Implementation structure

- Computes multipole moments from a set of points with user-specified
- spin populations. Syntax:
- Ilm=points2mult(xyz,mxyz,rho,L,method)
- The multipoles in question are described in Equation 32 of:
- xyz -coordinates as [x y z] with multiple rows,
- at which density rho is evaluated, in Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z], in
- Angstroms.
- L -array of ranks of spherical harmonics of the
- probability density
- rho -column of the densities at the points xyz
- method -'points' or 'grid'. If the spin density is supplied
