# kernel/optimcon/aux_mat.m

- Signature: `[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...`

## Purpose

Builds auxiliary matrices for the calculation of the directional derivatives of the trapezium product quadrature propagator: expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt) with respect to control coefficients in the evolution generators HL and HR on the left and the right edge of the interval. The de- rivatives are calculated using Eq 16 of Goodwin and Kuprov:

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...
cc_comm,dt,cL,cR,k,j)
```

## Parameters / inputs

- drifts -a cell array of two matrices containing drift
- generators at the left (first element) and the
- right (second element) edge of the interval
- controls -a cell array of K control generators
- cc_comm_idx -a KxK matrix of logicals indicating non-zero
- commutation of controls
- cc_comm -a KxK cell array control commutation relarions
- dt -interval duration, seconds
- cL -control generator coefficients at the left
- edge of the interval
- cR -control generator coefficients at the right
- edge of the interval
- k -the index of the generator inside controls
- array that the differentiation refers to
- j -(optional) the index of the 2nd generator
- inside controls array that the differentiation
- refers to. Required for 3x3 block auxiliary
- matrices

## Outputs

- auxm_l -auxilary matrix for the derivative of the
- interval propagator with respect to cL
- auxm_r -auxilary matrix for the derivative of the
- interval propagator with respect to cR

## Implementation structure

- Builds auxiliary matrices for the calculation of the directional
- derivatives of the trapezium product quadrature propagator:
- expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt)
- with respect to control coefficients in the evolution generators
- HL and HR on the left and the right edge of the interval. The de-
- rivatives are calculated using Eq 16 of Goodwin and Kuprov:
- [auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...
- cc_comm,dt,cL,cR,k,j)
- drifts -a cell array of two matrices containing drift
- generators at the left (first element) and the
- right (second element) edge of the interval
- controls -a cell array of K control generators
