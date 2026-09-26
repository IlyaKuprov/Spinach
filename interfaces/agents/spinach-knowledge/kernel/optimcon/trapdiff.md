# kernel/optimcon/trapdiff.m

- Signature: `[DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)`

## Purpose

Directional derivatives for the trapezium product quadrature publi- shed by Iserles and Norsett (see Corollary 3.3) in The derivatives are of the following propagator: expm(-i*((HL+HR)/2+i*dt*(sqrt(3)/12)*[HL,HR])*dt) with respect to the coefficients cL,cR in the evolution generators HL and HR on the left and the right side of the interval respecti- vely. Evolution generators HL and HR are split into the drift part H

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
[DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)
```

## Parameters / inputs

- Hd -a cell array of two matrices containing drift
- generators at the left (first element) and the
- right (second element) edge of the interval
- Hc -control operator or superoperator
- dt -interval duration, seconds
- cL -control operator coefficient at the
- left edge of the interval
- cR -control operator coefficient at the
- right edge of the interval

## Outputs

- DL -derivative of the interval propagator
- with respect to cL
- DR -derivative of the interval propagator
- with respect to cR

## Implementation structure

- Directional derivatives for the trapezium product quadrature publi-
- shed by Iserles and Norsett (see Corollary 3.3) in
- The derivatives are of the following propagator:
- expm(-i*((HL+HR)/2+i*dt*(sqrt(3)/12)*[HL,HR])*dt)
- with respect to the coefficients cL,cR in the evolution generators
- HL and HR on the left and the right side of the interval respecti-
- vely. Evolution generators HL and HR are split into the drift part
- Ho and the control part Hc, such that HL=Ho+cL*Hc and HR=Ho+cR*Hc
- on the left and the right edge of the interval.
- The derivatives are calculated using Eq 16 of Goodwin and Kuprov
- [DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)
- Hd -a cell array of two matrices containing drift
