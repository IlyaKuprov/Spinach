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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(l,m,theta,phi)`.
- Lines 29-30: Get Schmidt-normalized Legendres; implemented by `S=legendre(l,cos(theta),'sch')`.
- Lines 35-36: Make spherical harmonics; implemented by `if m==0`.
- Lines 42-43: Flip the sign if needed; implemented by `if (m>0)&&(mod(m,2)==1), Y=-Y; end`.

### Control flow inferred from the code

- Line 36: conditional branch on `m==0`.
- Line 43: conditional branch on `(m>0)&&(mod(m,2)==1), Y=-Y; end`.

### Key state/data transformations

- Lines 30: computes `S` using `S=legendre(l,cos(theta),'sch')`.
- Lines 37: computes `Y` using `Y=sqrt((2*l+1)/(4*pi))*S.*exp(1i*m*phi)`.

### Local helper functions

- Line 48: `grumble()` — `function grumble(l,m,theta,phi)`.
  - Representative operation: `if (~isnumeric(l))||(~isreal(l))||(~isscalar(l))||(mod(l,1)~=0)||(l<0)`.
  - Representative operation: `error('l must be a non-negative real integer.')`.

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
