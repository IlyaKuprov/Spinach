# kernel/plotting/efg_display.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/efg_display.m`
- Signature: `efg_display(props,atoms,scaling,conmatrix,options)`
- Total lines: 308

## Purpose

Electric field gradient tensors and their eigensystems. Two styles are implemented: A. Ellipsoids (symmetric tensors only): 1. A unit sphere in a Cartesian space is scaled by abs(Axx) in the x direction, abs(Ayy) in the y direction and abs(Azz) in the z direction, where Axx, Ayy, Azz are the eigenvalues of the CST ten- sor in units of ppm. 2. A set of axes is drawn inside the sphere with a red axis for a positive eig

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 71-72: Check consistency; implemented by `grumble(props,atoms,scaling,conmatrix)`.
- Lines 74-76: Set defaults; implemented by `if (~exist('options','var'))|| (~isfield(options,'style'))`.
- Lines 92-93: Set up graphics; implemented by `light('Position',[-2, 2,20])`.
- Lines 96-97: Draw the molecule; implemented by `molplot(props.std_geom,conmatrix)`.
- Lines 99-100: Sample a unit sphere; implemented by `k=5; npts=2^k-1`.
- Lines 108-109: Atoms by number or name; implemented by `if isnumeric(atoms)`.
- Lines 117-118: Loop over atoms; implemented by `for n=idx`.
- Lines 120-121: Get the current EFG or NQI tensor; implemented by `if isfield(props,'nqi')&&(~isempty(props.nqi{n}))`.
- Lines 129-130: Eliminate the isotroic part; implemented by `if options.kill_iso`.
- Lines 134-135: Pick the style; implemented by `switch options.style`.
- Lines 137-138: Ellipsoid plots; implemented by `case 'ellipsoids'`.
- Lines 140-141: Diagonalize the tensor; implemented by `[eigvecs,eigvals]=eig(efg,'vector')`.
- Lines 143-144: Make sure the tensor is symmetric; implemented by `if norm(eigvecs'*eigvecs-eye(3),2)>1e-3`.
- Lines 148-149: Scale vertex coordinates; implemented by `coords=[X*eigvals(1)*scaling`.
- Lines 153-154: Rotate vertex coordinates; implemented by `coords=eigvecs*coords`.
- Lines 156-157: Translate vertex coordinates; implemented by `coords=coords+props.std_geom(n,:)'`.
- Lines 159-160: Reshape for plotting; implemented by `Xp=reshape(coords(1,:),[npts npts])`.
- Lines 164-165: Simple grey colour; implemented by `Cp=ones([npts npts 3])/2`.

### Control flow inferred from the code

- Line 75: conditional branch on `(~exist('options','var'))||`.
- Line 79: conditional branch on `(~exist('options','var'))||`.
- Line 83: conditional branch on `(~exist('options','var'))||`.
- Line 87: conditional branch on `(~exist('options','var'))||`.
- Line 109: conditional branch on `isnumeric(atoms)`.
- Line 118: `for` loop over `n=idx`.
- Line 121: conditional branch on `isfield(props,'nqi')&&(~isempty(props.nqi{n}))`.
- Line 130: conditional branch on `options.kill_iso`.
- Line 135: dispatches on `options.style`; cases `'ellipsoids'`, `'harmonics'`.
- Line 144: conditional branch on `norm(eigvecs'*eigvecs-eye(3),2)>1e-3`.
- Line 178: conditional branch on `eigvals(1)<0, col_a='b-'; end`.
- Line 179: conditional branch on `eigvals(2)<0, col_b='b-'; end`.
- Line 180: conditional branch on `eigvals(3)<0, col_c='b-'; end`.
- Line 240: conditional branch on `options.numbers`.

### Key state/data transformations

- Lines 77: computes `options.style` using `options.style='harmonics'`.
- Lines 81: computes `options.kill_iso` using `options.kill_iso=false()`.
- Lines 85: computes `options.numbers` using `options.numbers=false()`.
- Lines 89: computes `options.symbols` using `options.symbols=true()`.
- Lines 100: computes `k` using `k=5; npts=2^k-1`.
- Lines 101: computes `theta` using `theta=linspace(0,pi,npts)`.
- Lines 102: computes `phi` using `phi=linspace(0,2*pi,npts)`.
- Lines 103: computes `[T,P]` using `[T,P]=ndgrid(theta,phi)`.
- Lines 104: computes `X` using `X=sin(T).*cos(P); X=X(:)'`.
- Lines 105: computes `Y` using `Y=sin(T).*sin(P); Y=Y(:)'`.
- Lines 106: computes `Z` using `Z=cos(T); Z=Z(:)'`.
- Lines 110: computes `idx` using `idx=atoms(:)'`.
- Lines 131: computes `efg` using `efg=efg-eye(3)*trace(efg)/3`.
- Lines 141: computes `[eigvecs,eigvals]` using `[eigvecs,eigvals]=eig(efg,'vector')`.
- Lines 145: computes `error('significant antisymmetry found, use options.style` using `error('significant antisymmetry found, use options.style=''harmonics''')`.
- Lines 149: computes `coords` using `coords=[X*eigvals(1)*scaling`.
- Lines 160: computes `Xp` using `Xp=reshape(coords(1,:),[npts npts])`.
- Lines 161: computes `Yp` using `Yp=reshape(coords(2,:),[npts npts])`.

### Local helper functions

- Line 269: `grumble()` — `function grumble(props,atoms,scaling,conmatrix)`.
  - Representative operation: `if ~isfield(props,'std_geom')`.
  - Representative operation: `error('props structure does not contain std_geom field.')`.

## Syntax

```matlab
efg_display(props,atoms,scaling,conmatrix,options)
```

## Parameters / inputs

- props -output of c2spinach or gparse
- atoms -a cell array of element symbols
- or a vector of integers, indica-
- ting the atoms for which EFG
- tensors should be visualised,
- e.g. {'N','O'} or [1 2 5]
- scaling -a factor to scale the tensors
- by for visualisation
- conmatrix -binary connectivity matrix, 1
- if a pair of atoms should be
- connected by a bond. If an em-
- pty vector is supplied, 1.6
- Angstrom cutoff distance is used
- options.style -'ellipsoids' or 'harmonics'
- options.kill_iso -set to true() to eliminate the
- isotropic parts of tensors be-
- fore plotting
- options.numbers -set to true() to display atom
- numbers
- options.symbols -set to false() to not display
- atom symbols

## Implementation structure

- Electric field gradient tensors and their eigensystems. Two
- styles are implemented:
- A. Ellipsoids (symmetric tensors only):
- 1. A unit sphere in a Cartesian space is scaled by
- abs(Axx) in the x direction, abs(Ayy) in the y
- direction and abs(Azz) in the z direction, where
- Axx, Ayy, Azz are the eigenvalues of the CST ten-
- sor in units of ppm.
- 2. A set of axes is drawn inside the sphere with a
- red axis for a positive eigenvalue, and a blue
- axis for a negative one.
- 3. The sphere is translated to the point of corres-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `exist()`, `isfield()`, `false()`, `true()`, `light()`, `molplot()`, `atoms()`, `iscell()`, `ismember()`, `num2str()`, `eigvals()`, `coords()`, `eigvecs()`, `plot3()`, `vector_a()`.
