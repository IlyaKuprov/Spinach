# kernel/optimcon/trapdiff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/trapdiff.m`
- Signature: `[DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)`
- Total lines: 97

## Purpose

Directional derivatives for the trapezium product quadrature publi- shed by Iserles and Norsett (see Corollary 3.3) in The derivatives are of the following propagator: expm(-i*((HL+HR)/2+i*dt*(sqrt(3)/12)*[HL,HR])*dt) with respect to the coefficients cL,cR in the evolution generators HL and HR on the left and the right side of the interval respecti- vely. Evolution generators HL and HR are split into the drift part H

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 55-56: Check consistency; implemented by `grumble(Hd,Hc,dt,cL,cR)`.
- Lines 58-59: Precompute directions; implemented by `H_dir_L=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hc*Hd{2}-Hd{2}*Hc)`.
- Lines 62-63: Call directional derivative function; implemented by `D=dirdiff(spin_system,(Hd{1}+Hd{2})/2+cL*H_dir_L+cR*H_dir_R,H_dir_L,dt,2); DL=D{2}`.

### Key state/data transformations

- Lines 59: computes `H_dir_L` using `H_dir_L=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hc*Hd{2}-Hd{2}*Hc)`.
- Lines 60: computes `H_dir_R` using `H_dir_R=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hd{1}*Hc-Hc*Hd{1})`.
- Lines 63: computes `D` using `D=dirdiff(spin_system,(Hd{1}+Hd{2})/2+cL*H_dir_L+cR*H_dir_R,H_dir_L,dt,2); DL=D{2}`.

### Local helper functions

- Line 69: `grumble()` — `function grumble(Hd,Hc,dt,cL,cR)`.
  - Representative operation: `if ~iscell(Hd), error('Hd must be a cell array of matrices.'); end`.
  - Representative operation: `if numel(Hd)~=2`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dirdiff()`, `iscell()`, `isscalar()`.
