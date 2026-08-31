# kernel/conventions/transforms/cart2mode.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/cart2mode.m`
- Signature: `mode_derivs=cart2mode(cart_derivs,eigvecs,masses,frqs)`
- Total lines: 131

## Purpose

Converts Cartesian derivatives of spin Hamiltonian parameters, as produced by electronic structure theory packages, into the derivatives with respect to dimensionless mode coordinates that the bosonic mode specification interface of create.m expects in inter.modes.coupling_mod and inter.modes.zeeman_mod fields. The dimensionless coordinate of each mode is (a+a')/sqrt(2), and the Cartesian displacement that it produce

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 58-59: Check consistency; implemented by `grumble(cart_derivs,eigvecs,masses,frqs)`.
- Lines 61-62: Expand atomic masses over Cartesian degrees of freedom; implemented by `dof_masses=kron(masses,[1;1;1])*1.66053906892e-27`.
- Lines 64-65: Zero-point displacement scale vectors in Angstrom; implemented by `scales=1e10*eigvecs.*sqrt((6.62607015e-34/(2*pi))./(dof_masses*(2*pi*frqs)))`.
- Lines 67-68: Contract the derivatives with the scale vectors; implemented by `if ndims(cart_derivs)==3`.

### Control flow inferred from the code

- Line 68: conditional branch on `ndims(cart_derivs)==3`.

### Key state/data transformations

- Lines 62: computes `dof_masses` using `dof_masses=kron(masses,[1;1;1])*1.66053906892e-27`.
- Lines 65: computes `scales` using `scales=1e10*eigvecs.*sqrt((6.62607015e-34/(2*pi))./(dof_masses*(2*pi*frqs)))`.
- Lines 69: computes `mode_derivs` using `mode_derivs=sum(cart_derivs.*reshape(scales,1,1,[]),3)`.

### Local helper functions

- Line 78: `grumble()` — `function grumble(cart_derivs,eigvecs,masses,frqs)`.
  - Representative operation: `if (~isnumeric(cart_derivs))||(~isreal(cart_derivs))|| any(~isfinite(cart_derivs),'all')`.
  - Representative operation: `any(~isfinite(cart_derivs),'all')`.

## Parameters / inputs

- cart_derivs -first derivatives of an interaction with
- respect to Cartesian displacements, in Hz
- per Angstrom, an array of dimension
- [d1 d2 3N] where N is the number of atoms
- and the Cartesian degrees of freedom are
- ordered [x1 y1 z1 x2 y2 z2 ...]; or second
- derivatives in Hz per Angstrom squared, an
- array of dimension [d1 d2 3N 3N]
- eigvecs -orthonormal mass-weighted normal mode
- eigenvector, a [3N 1] column vector for
- the first order case; two such vectors
- as a [3N 2] array for the second order
- case
- masses -atomic masses in unified atomic mass
- units, an [N 1] column vector
- frqs -mode frequency in Hz, a positive scalar
- for the first order case; a [1 2] vector
- for the second order case

## Outputs

- mode_derivs -derivatives with respect to dimensionless
- mode coordinates, in Hz, a [d1 d2] array
- to be placed into the corresponding cell
- of inter.modes.coupling_mod (d1=3, d2=3)
- or inter.modes.zeeman_mod (d1=1, d2=3)
- Note: raw Taylor derivatives are returned; the 1/2 factors of
- the Taylor expansion are applied by Spinach internally.
- Derivative data in wavenumbers or meV should be conver-
- ted into Hz with icm2hz.m or mev2hz.m beforehand. Zero
- and negative frequency modes are rejected because their
- zero-point scaling is undefined.

## Implementation structure

- Converts Cartesian derivatives of spin Hamiltonian parameters,
- as produced by electronic structure theory packages, into the
- derivatives with respect to dimensionless mode coordinates that
- the bosonic mode specification interface of create.m expects in
- inter.modes.coupling_mod and inter.modes.zeeman_mod fields. The
- dimensionless coordinate of each mode is (a+a')/sqrt(2), and the
- Cartesian displacement that it produces on a degree of freedom
- with mass m is scaled by sqrt(hbar/(m*omega)), where omega is
- the angular frequency of the mode. Syntax:
- mode_derivs=cart2mode(cart_derivs,eigvecs,masses,frqs)
- cart_derivs -first derivatives of an interaction with
- respect to Cartesian displacements, in Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ndims()`, `scales()`, `any()`, `three()`, `four()`, `iscolumn()`, `isscalar()`, `isequal()`.
