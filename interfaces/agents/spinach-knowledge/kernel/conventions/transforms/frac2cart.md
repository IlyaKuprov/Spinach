# kernel/conventions/transforms/frac2cart.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/frac2cart.m`
- Signature: `[XYZ,va,vb,vc]=frac2cart(a,b,c,alp,bet,gam,ABC)`
- Total lines: 67

## Purpose

Converts fractional crystallographic coordinates to Cartesian coordinates. Syntax: [XYZ,va,vb,vc]=frac2cart(a,b,c,alpha,beta,gamma,ABC)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(a,b,c,alp,bet,gam,ABC)`.
- Lines 31-32: Compute the transformation matrix; implemented by `v=a*b*c*sqrt(1-cosd(alp)^2-cosd(bet)^2-cosd(gam)^2+2*cosd(alp)*cosd(bet)*cosd(gam))`.
- Lines 37-38: Apply the transformation matrix; implemented by `XYZ=(T*ABC')'`.
- Lines 40-41: Get the primitive vectors; implemented by `va=T(:,1); vb=T(:,2); vc=T(:,3)`.

### Key state/data transformations

- Lines 32: computes `v` using `v=a*b*c*sqrt(1-cosd(alp)^2-cosd(bet)^2-cosd(gam)^2+2*cosd(alp)*cosd(bet)*cosd(gam))`.
- Lines 33: computes `T` using `T=[a b*cosd(gam) c*cosd(bet)`.
- Lines 38: computes `XYZ` using `XYZ=(T*ABC')'`.
- Lines 41: computes `va` using `va=T(:,1); vb=T(:,2); vc=T(:,3)`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(a,b,c,alp,bet,gam,ABC)`.
  - Representative operation: `if (~isnumeric(a))||(~isscalar(a))||(~isreal(a))||(a<=0)|| (~isnumeric(b))||(~isscalar(b))||(~isreal(b))||(b<=0)|| (~isnumeric(c))||(~isscalar(c))||(~isreal(c))||(c<=0)`.
  - Representative operation: `(~isnumeric(b))||(~isscalar(b))||(~isreal(b))||(b<=0)|| (~isnumeric(c))||(~isscalar(c))||(~isreal(c))||(c<=0)`.

## Parameters / inputs

- a,b,c -three unit cell dimensions
- alp,bet,gam -three unit cell angles, degrees
- ABC -fractional atomic coordinates as
- Nx3 array of numbers

## Outputs

- XYZ -Cartesian atomic coordinates as
- Nx3 array of numbers
- va, vb, vc -primitive lattice vectors

## Implementation structure

- Converts fractional crystallographic coordinates to Cartesian
- coordinates. Syntax:
- [XYZ,va,vb,vc]=frac2cart(a,b,c,alpha,beta,gamma,ABC)
- a,b,c -three unit cell dimensions
- alp,bet,gam -three unit cell angles, degrees
- ABC -fractional atomic coordinates as
- Nx3 array of numbers
- XYZ -Cartesian atomic coordinates as
- va, vb, vc -primitive lattice vectors
- Check consistency
- Compute the transformation matrix
- Apply the transformation matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cosd()`, `sind()`, `isscalar()`.
