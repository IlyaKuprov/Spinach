# kernel/derivatives/fftdiff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fftdiff.m`
- Signature: `kern=fftdiff(order,npoints,dx)`
- Total lines: 74

## Purpose

Spectral differentiation kernel. Syntax: kern=fftdiff(order,npoints,dx)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(order,npoints,dx)`.
- Lines 32-33: Adapt to the point count; implemented by `if mod(npoints,2)==1`.
- Lines 35-36: Kernel for odd point counts; implemented by `kern=ifftshift((2i*pi*((1-npoints)/2:((npoints)/2))/(npoints*dx)).^order)`.
- Lines 40-41: Kernel for even point counts; implemented by `kern=ifftshift((2i*pi*(((-npoints)/2):((npoints-1)/2))/(npoints*dx)).^order)`.
- Lines 45-46: Complain and bomb out; implemented by `error('npoints parameter must be an integer.')`.

### Control flow inferred from the code

- Line 33: conditional branch on `mod(npoints,2)==1`.

### Key state/data transformations

- Lines 36: computes `kern` using `kern=ifftshift((2i*pi*((1-npoints)/2:((npoints)/2))/(npoints*dx)).^order)`.

### Local helper functions

- Line 53: `grumble()` — `function grumble(order,npoints,dx)`.
  - Representative operation: `if (~isnumeric(order))||(~isreal(order))||(numel(order)~=1)|| (order<1)||(mod(order,1)~=0)`.
  - Representative operation: `(order<1)||(mod(order,1)~=0)`.

## Parameters / inputs

- order -order of the derivative
- npoints -number of points in the grid
- dx -grid step length
- Output:
- kern -a vector that is to be used for accurate
- numerical differentiation of periodic real
- signals in the following way:
- derivative=real(ifft(fft(signal).*kern));
- Note: periodic boundary conditions.

## Implementation structure

- Spectral differentiation kernel. Syntax:
- kern=fftdiff(order,npoints,dx)
- order -order of the derivative
- npoints -number of points in the grid
- dx -grid step length
- Output:
- kern -a vector that is to be used for accurate
- numerical differentiation of periodic real
- signals in the following way:
- derivative=real(ifft(fft(signal).*kern));
- Note: periodic boundary conditions.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ifftshift()`.
