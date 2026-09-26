# kernel/derivatives/fftdiff.m

- Signature: `kern=fftdiff(order,npoints,dx)`

## Purpose

Spectral differentiation kernel. Syntax: kern=fftdiff(order,npoints,dx)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

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
