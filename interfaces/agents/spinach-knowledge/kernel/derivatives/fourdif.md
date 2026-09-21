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
