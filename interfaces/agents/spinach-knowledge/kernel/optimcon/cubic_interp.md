# kernel/optimcon/cubic_interp.m

- Signature: `[alpha,fx]=cubic_interp(end_a,end_b,alpha_a,alpha_b,...`

## Purpose

Finds the extremum of a cubic interpolant built from function values and directional derivatives at two points and returns the best point inside the interpolation interval. Syntax: [alpha,fx]=cubic_interp(end_A,end_B,alpha_A,alpha_B,... f_A,dir_deriv_A,f_B,dir_deriv_B)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

## Parameters / inputs

- end_a -first interpolation boundary in alpha space
- end_b -second interpolation boundary in alpha space
- alpha_a -first interpolation anchor point
- alpha_b -second interpolation anchor point
- f_a -function value at alpha_a
- dir_der_a -directional derivative at alpha_a
- f_b -function value at alpha_b
- dir_der_b -directional derivative at alpha_b

## Outputs

- alpha -selected maximiser of the cubic model
- fx -cubic model value at alpha

## Implementation structure

- Finds the extremum of a cubic interpolant built from function
- values and directional derivatives at two points and returns
- the best point inside the interpolation interval. Syntax:
- [alpha,fx]=cubic_interp(end_A,end_B,alpha_A,alpha_B,...
- f_A,dir_deriv_A,f_B,dir_deriv_B)
- end_a -first interpolation boundary in alpha space
- end_b -second interpolation boundary in alpha space
- alpha_a -first interpolation anchor point
- alpha_b -second interpolation anchor point
- f_a -function value at alpha_a
- dir_der_a -directional derivative at alpha_a
- f_b -function value at alpha_b
