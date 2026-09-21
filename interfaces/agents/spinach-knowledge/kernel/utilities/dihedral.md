# kernel/utilities/dihedral.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dihedral.m`
- Signature: `phi=dihedral(A,B,C,D)`
- Total lines: 65

## Purpose

Computes the dihedral angle between vectors specified by the four sets of atomic coordinates. The atoms are assu- med to be bonded as A-B-C-D. Syntax: phi=dihedral(A,B,C,D)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `atan2()`, `dot()`, `cross()`, `isrow()`.
