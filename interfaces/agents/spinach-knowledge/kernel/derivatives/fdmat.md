# kernel/derivatives/fdmat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fdmat.m`
- Signature: `D=fdmat(dim,nstenc,order,boundary)`
- Total lines: 100

## Purpose

Returns arbitrary-order central finite-difference differentiation matrices (sparse) with unit grid point spacing. Syntax: D=fdmat(dim,nstenc,order,boundary)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Set the default; implemented by `if ~exist('boundary','var'), boundary='pbc'; end`.
- Lines 34-35: Check consistency; implemented by `grumble(dim,nstenc,order,boundary)`.
- Lines 37-38: Preallocate the answer; implemented by `D=spalloc(dim,dim,dim*nstenc)`.
- Lines 40-41: Decide the filling procedure; implemented by `switch boundary`.
- Lines 45-46: Edges filled with sided schemes; implemented by `for n=1:(nstenc-1)/2`.
- Lines 51-52: Middle filled with centered schemes; implemented by `stencil=((1-nstenc)/2):((nstenc-1)/2)`.
- Lines 60-61: Wraparound fill with centered schemes; implemented by `stencil=((1-nstenc)/2):((nstenc-1)/2)`.
- Lines 69-70: Complain and bomb out; implemented by `error('unknown boundary type.')`.

### Control flow inferred from the code

- Line 32: conditional branch on `~exist('boundary','var'), boundary='pbc'; end`.
- Line 41: dispatches on `boundary`; cases `'wall'`, `'pbc'`.
- Line 46: `for` loop over `n=1:(nstenc-1)/2`.
- Line 54: `for` loop over `n=((nstenc-1)/2+1):(dim-(nstenc-1)/2)`.
- Line 63: `for` loop over `n=1:dim`.

### Key state/data transformations

- Lines 38: computes `D` using `D=spalloc(dim,dim,dim*nstenc)`.
- Lines 47: computes `w` using `w=fdweights(n,1:nstenc,order); D(n,1:nstenc)=w(end,:)`.
- Lines 48: computes `D(end-n+1,(end-nstenc+1):end)` using `D(end-n+1,(end-nstenc+1):end)=((-1)^order)*w(end,end:-1:1)`.
- Lines 52: computes `stencil` using `stencil=((1-nstenc)/2):((nstenc-1)/2)`.
- Lines 55: computes `D(n,stencil+n)` using `D(n,stencil+n)=w(end,:)`.
- Lines 64: computes `D(n,mod(stencil+n-1,dim)+1)` using `D(n,mod(stencil+n-1,dim)+1)=w(end,:)`.

### Local helper functions

- Line 77: `grumble()` — `function grumble(dim,nstenc,order,boundary)`.
  - Representative operation: `if (dim<1)||(nstenc<1)||(order<1)||(mod(dim,1)~=0)|| (mod(nstenc,1)~=0)||(mod(order,1)~=0)`.
  - Representative operation: `(mod(nstenc,1)~=0)||(mod(order,1)~=0)`.

## Parameters / inputs

- dim -dimension of the column vector to be
- differentiated
- nstenc -number of points in the finite diffe-
- rence stencil
- order -order of the derivative required
- boundary -'wall' fills the edges with sided
- finite difference schemes, 'pbc'
- assumes periodic boundaries. The
- default is 'pbc'.

## Outputs

- D -finite difference differentiation matrix

## Implementation structure

- Returns arbitrary-order central finite-difference differentiation
- matrices (sparse) with unit grid point spacing. Syntax:
- D=fdmat(dim,nstenc,order,boundary)
- dim -dimension of the column vector to be
- differentiated
- nstenc -number of points in the finite diffe-
- rence stencil
- order -order of the derivative required
- boundary -'wall' fills the edges with sided
- finite difference schemes, 'pbc'
- assumes periodic boundaries. The
- default is 'pbc'.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `spalloc()`, `fdweights()`, `ischar()`.
