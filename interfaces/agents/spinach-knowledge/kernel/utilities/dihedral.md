# kernel/utilities/dihedral.m

- Signature: `phi=dihedral(A,B,C,D)`

## Purpose

Computes the dihedral angle between vectors specified by the four sets of atomic coordinates. The atoms are assu- med to be bonded as A-B-C-D. Syntax: phi=dihedral(A,B,C,D)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- A -row vector of cartesian coordinates
- for atom A
- B -row vector of cartesian coordinates
- for atom B
- C -row vector of cartesian coordinates
- for atom C
- D -row vector of cartesian coordinates
- for atom D

## Outputs

- phi -dihedral angle, degrees

## Implementation structure

- Computes the dihedral angle between vectors specified by
- the four sets of atomic coordinates. The atoms are assu-
- med to be bonded as A-B-C-D. Syntax:
- phi=dihedral(A,B,C,D)
- A - row vector of cartesian coordinates
- for atom A
- B - row vector of cartesian coordinates
- for atom B
- C - row vector of cartesian coordinates
- for atom C
- D - row vector of cartesian coordinates
- for atom D
