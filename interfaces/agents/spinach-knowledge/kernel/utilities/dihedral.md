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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(A,B,C,D)`.
- Lines 34-35: Do the math; implemented by `b1=(B-A)/norm(B-A,2)`.

### Key state/data transformations

- Lines 35: computes `b1` using `b1=(B-A)/norm(B-A,2)`.
- Lines 36: computes `b2` using `b2=(C-B)/norm(C-B,2)`.
- Lines 37: computes `b3` using `b3=(D-C)/norm(D-C,2)`.
- Lines 38-39: computes `phi` using `phi=180*atan2(dot(norm(b2,2)*b1,cross(b2,b3)), dot(cross(b1,b2),cross(b2,b3)))/pi`.

### Local helper functions

- Line 44: `grumble()` — `function grumble(A,B,C,D)`. Why is the plight of the poor not a matter of more sustained public discussion and more decisive government policy? [...] For many Ameri-
  - Representative operation: `if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))||(~isnumeric(D))|| (~isreal(A))||(~isreal(B))||(~isreal(C))||(~isreal(D))|| (numel(A)~=3)||(numel(B)~=3)||(numel(C)~…`.
  - Representative operation: `(~isreal(A))||(~isreal(B))||(~isreal(C))||(~isreal(D))|| (numel(A)~=3)||(numel(B)~=3)||(numel(C)~=3)||(numel(D)~=3)|| (~isrow(A))||(~isrow(B))||(~isrow(C))||(~isrow(D))`.

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
