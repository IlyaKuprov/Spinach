# kernel/utilities/conmat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/conmat.m`
- Signature: `conmatrix=conmat(xyz,r0)`
- Total lines: 106

## Purpose

Molecular connectivity matrix calculator with N*log(N) asymptotic complexity scaling with respect to the num- ber or atoms. Syntax: conmatrix=conmat(xyz,r0)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(xyz,r0)`.
- Lines 31-32: Sort X coordinates; implemented by `[x_sorted,x_index]=sort(xyz(:,1))`.
- Lines 34-35: Scan X coordinates; implemented by `A=false(size(xyz,1),size(xyz,1))`.
- Lines 47-48: Sort Y coordinates; implemented by `[y_sorted,y_index]=sort(xyz(:,2))`.
- Lines 50-51: Scan Y coordinates; implemented by `B=false(size(xyz,1),size(xyz,1))`.
- Lines 63-64: Sort Z coordinates; implemented by `[z_sorted,z_index]=sort(xyz(:,3))`.
- Lines 66-67: Scan Y coordinates; implemented by `C=false(size(xyz,1),size(xyz,1))`.
- Lines 79-80: Compile the dirty matrix; implemented by `conmatrix=A&B&C`.
- Lines 82-83: Clean up the matrix; implemented by `[row,col]=find(conmatrix)`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=1:size(xyz,1)`.
- Line 37: `for` loop over `k=(n+1):size(xyz,1)`.
- Line 38: conditional branch on `abs(x_sorted(n)-x_sorted(k))<r0`.
- Line 52: `for` loop over `n=1:size(xyz,1)`.
- Line 53: `for` loop over `k=(n+1):size(xyz,1)`.
- Line 54: conditional branch on `abs(y_sorted(n)-y_sorted(k))<r0`.
- Line 68: `for` loop over `n=1:size(xyz,1)`.
- Line 69: `for` loop over `k=(n+1):size(xyz,1)`.
- Line 70: conditional branch on `abs(z_sorted(n)-z_sorted(k))<r0`.
- Line 84: `for` loop over `n=1:numel(row)`.
- Line 85: conditional branch on `norm(xyz(row(n),:)-xyz(col(n),:),2)>r0`.

### Key state/data transformations

- Lines 32: computes `[x_sorted,x_index]` using `[x_sorted,x_index]=sort(xyz(:,1))`.
- Lines 35: computes `A` using `A=false(size(xyz,1),size(xyz,1))`.
- Lines 39: computes `A(x_index(n),x_index(k))` using `A(x_index(n),x_index(k))=1`.
- Lines 40: computes `A(x_index(k),x_index(n))` using `A(x_index(k),x_index(n))=1`.
- Lines 48: computes `[y_sorted,y_index]` using `[y_sorted,y_index]=sort(xyz(:,2))`.
- Lines 51: computes `B` using `B=false(size(xyz,1),size(xyz,1))`.
- Lines 55: computes `B(y_index(n),y_index(k))` using `B(y_index(n),y_index(k))=1`.
- Lines 56: computes `B(y_index(k),y_index(n))` using `B(y_index(k),y_index(n))=1`.
- Lines 64: computes `[z_sorted,z_index]` using `[z_sorted,z_index]=sort(xyz(:,3))`.
- Lines 67: computes `C` using `C=false(size(xyz,1),size(xyz,1))`.
- Lines 71: computes `C(z_index(n),z_index(k))` using `C(z_index(n),z_index(k))=1`.
- Lines 72: computes `C(z_index(k),z_index(n))` using `C(z_index(k),z_index(n))=1`.
- Lines 80: computes `conmatrix` using `conmatrix=A&B&C`.
- Lines 83: computes `[row,col]` using `[row,col]=find(conmatrix)`.
- Lines 86: computes `conmatrix(row(n),col(n))` using `conmatrix(row(n),col(n))=0`.

### Local helper functions

- Line 94: `grumble()` — `function grumble(xyz,r0)`. It's always good to be underestimated. Donald Trump
  - Representative operation: `if (~isnumeric(xyz))||(size(xyz,2)~=3)||(~isreal(xyz))`.
  - Representative operation: `error('xyz parameter should be a real matrix with three columns.')`.

## Parameters / inputs

- xyz -an array with N rows and three columns,
- giving the Cartesian coordinates of
- each particle
- r0 -the distance below which the particles
- are to be considered "connected"
- Output:
- conmatrix -a sparse logical matrix containing 1
- at the positions corresponding to the
- connected particles

## Implementation structure

- Molecular connectivity matrix calculator with N*log(N)
- asymptotic complexity scaling with respect to the num-
- ber or atoms. Syntax:
- conmatrix=conmat(xyz,r0)
- xyz -an array with N rows and three columns,
- giving the Cartesian coordinates of
- each particle
- r0 -the distance below which the particles
- are to be considered "connected"
- Output:
- conmatrix -a sparse logical matrix containing 1
- at the positions corresponding to the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xyz()`, `false()`, `x_sorted()`, `x_index()`, `y_sorted()`, `y_index()`, `z_sorted()`, `z_index()`, `row()`, `col()`, `conmatrix()`.
