# kernel/utilities/spher_harmon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/spher_harmon.m`
- Signature: `Y=spher_harmon(l,m,theta,phi)`
- Total lines: 67

## Purpose

Spherical harmonics. Syntax: Y=spher_harmon(l,m,theta,phi)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- l -L quantum number
- m -M quantum number
- theta -an array of theta angles in radians
- phi -an array of phi angles in radians

## Outputs

- Y -an array of spherical harmonics
- evaluated at the angles specified

## Implementation structure

- Spherical harmonics. Syntax:
- Y=spher_harmon(l,m,theta,phi)
- l -L quantum number
- m -M quantum number
- theta -an array of theta angles in radians
- phi -an array of phi angles in radians
- Y -an array of spherical harmonics
- evaluated at the angles specified
- Check consistency
- Get Schmidt-normalized Legendres
- Make spherical harmonics
- Flip the sign if needed

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `legendre()`, `squeeze()`, `isscalar()`.
