# kernel/derivatives/fourdif.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fourdif.m`
- Signature: `[x,DM]=fourdif(N,m)`
- Total lines: 98

## Purpose

The function [x,DM] = fourdif(N,m) computes the m'th derivative Fourier spectral differentiation matrix on grid with N equispa- ced points in [0,2pi). Syntax: [x,DM]=fourdif(N,m)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(N,m)`.
- Lines 37-38: Run the blob; implemented by `x=2*pi*(0:N-1)'/N; % grid points`.
- Lines 40-41: Differentiation matrix; implemented by `if (nargout>1)`.

### Control flow inferred from the code

- Line 41: conditional branch on `(nargout>1)`.
- Line 46: conditional branch on `m==0`.
- Line 50: conditional branch on `rem(N,2)==0`.
- Line 59: conditional branch on `rem(N,2)==0`.
- Line 72: conditional branch on `rem(m,2)==0`.

### Key state/data transformations

- Lines 38: computes `x` using `x=2*pi*(0:N-1)'/N; % grid points`.
- Lines 42: computes `h` using `h=2*pi/N`.
- Lines 43: computes `kk` using `kk=(1:N-1)'`.
- Lines 44: computes `n1` using `n1=floor((N-1)/2)`.
- Lines 45: computes `n2` using `n2=ceil((N-1)/2)`.
- Lines 47: computes `col1` using `col1=[1; zeros(N-1,1)]`.
- Lines 48: computes `row1` using `row1=col1`.
- Lines 51: computes `topc` using `topc=cot((1:n2)'*h/2)`.
- Lines 68: computes `N1` using `N1=floor((N-1)/2)`.
- Lines 69: computes `N2` using `N2=(-N/2)*rem(m+1,2)*ones(rem(N+1,2))`.
- Lines 70: computes `mwave` using `mwave=1i*[(0:N1) N2 (-N1:-1)]`.
- Lines 79: computes `DM` using `DM=toeplitz(col1,row1)`.

### Local helper functions

- Line 84: `grumble()` — `function grumble(N,m)`. When you first kissed me, the desire to kiss you back has only very narrowly overpowered the desire to slap
  - Representative operation: `if (~isnumeric(N))||(~isreal(N))||(numel(N)~=1)||(N<1)||(mod(N,1)~=0)`.
  - Representative operation: `error('N must be a positive real integer')`.

## Parameters / inputs

- N -dimension of differentiation matrix
- m -derivative order

## Outputs

- x -equispaced points 0, 2*pi/N, 4*pi/N, ... , (N-1)*2*pi/N
- DM -m-th order differentiation matrix
- Explicit formulae are used to compute the matrices for m=1 and
- m=2. A discrete Fourier approach is employed for m>2. The prog-
- ram computes the first column and first row and then uses the
- toeplitz() function to create the matrix.
- For m=1 and 2 the code implements a "flipping trick" to improve
- accuracy as suggested in http://dx.doi.org/10.1137/0916073
- S.C. Reddy
- J.A.C. Weideman

## Implementation structure

- The function [x,DM] = fourdif(N,m) computes the m'th derivative
- Fourier spectral differentiation matrix on grid with N equispa-
- ced points in [0,2pi). Syntax:
- [x,DM]=fourdif(N,m)
- N -dimension of differentiation matrix
- m -derivative order
- x -equispaced points 0, 2*pi/N, 4*pi/N, ... , (N-1)*2*pi/N
- DM -m-th order differentiation matrix
- Explicit formulae are used to compute the matrices for m=1 and
- m=2. A discrete Fourier approach is employed for m>2. The prog-
- ram computes the first column and first row and then uses the
- toeplitz() function to create the matrix.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cot()`, `flipud()`, `topc()`, `csc()`, `col1()`, `toeplitz()`.
