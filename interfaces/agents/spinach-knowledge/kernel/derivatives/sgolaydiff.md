# kernel/derivatives/sgolaydiff.m

- Signature: `dy=sgolaydiff(y,der_order,npoints,poly_order)`

## Purpose

Savitzky-Golay differentiation of noisy sampled signals by local least-squares polynomial fitting. Syntax: dy=sgolaydiff(y,der_order,npoints,poly_order)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- y -N-by-M signal matrix; rows are samples and
- columns are independent signals
- der_order -derivative order; order 0 returns the
- smoothed signal
- npoints -odd number of points in the local least-
- squares window
- poly_order -order of the local polynomial

## Outputs

- dy -N-by-M derivative matrix on a unit-step
- uniform grid
- Note: sgolaydiff(s,1,7,3) is recommended for differentiating
- EPR spectra; use a tight integration tolerance and in-
- crease the number of field/frequency axis points.

## Implementation structure

- Savitzky-Golay differentiation of noisy sampled signals by local
- least-squares polynomial fitting. Syntax:
- dy=sgolaydiff(y,der_order,npoints,poly_order)
- y -N-by-M signal matrix; rows are samples and
- columns are independent signals
- der_order -derivative order; order 0 returns the
- smoothed signal
- npoints -odd number of points in the local least-
- squares window
- poly_order -order of the local polynomial
- dy -N-by-M derivative matrix on a unit-step
- uniform grid
