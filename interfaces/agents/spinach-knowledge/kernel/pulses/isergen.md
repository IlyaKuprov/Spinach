# kernel/pulses/isergen.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/isergen.m`
- Signature: `H=isergen(HL,HM,HR,dt)`
- Total lines: 72

## Purpose

2nd and 4th order Iserles product quadrature generators for one time propagation step in the case of state-inde- pendent Hamiltonian. Syntax: H=isergen(HL,HM,HR,dt)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(HL,HM,HR,dt)`.
- Lines 34-35: Decide the quadrature; implemented by `if isempty(HM)`.
- Lines 37-38: Second order product quadrature; implemented by `H=(HL+HR)/2+(1i*dt/6)*(HL*HR-HR*HL)`.
- Lines 42-43: Fourth order product quadrature; implemented by `H=(HL+4*HM+HR)/6+(1i*dt/12)*(HL*HR-HR*HL)`.

### Control flow inferred from the code

- Line 35: conditional branch on `isempty(HM)`.

### Key state/data transformations

- Lines 38: computes `H` using `H=(HL+HR)/2+(1i*dt/6)*(HL*HR-HR*HL)`.

### Local helper functions

- Line 50: `grumble()` — `function grumble(HL,HM,HR,dt)`.
  - Representative operation: `if (~isnumeric(HL))||(size(HL,1)~=size(HL,2))`.
  - Representative operation: `error('HL must be a square matrix')`.

## Parameters / inputs

- HL -Hamiltonian at the left edge of the interval
- HM -[optional] Hamiltonian at the interval mid-
- point; if this is empty, second order quad-
- rature is used.
- HR -Hamiltonian at the right edge of the interval
- dt -interval duration, seconds

## Outputs

- H -effective evolution generator, to be used
- as exp(-1i*H*dt)

## Implementation structure

- 2nd and 4th order Iserles product quadrature generators
- for one time propagation step in the case of state-inde-
- pendent Hamiltonian. Syntax:
- H=isergen(HL,HM,HR,dt)
- HL -Hamiltonian at the left edge of the interval
- HM -[optional] Hamiltonian at the interval mid-
- point; if this is empty, second order quad-
- rature is used.
- HR -Hamiltonian at the right edge of the interval
- dt -interval duration, seconds
- H -effective evolution generator, to be used
- as exp(-1i*H*dt)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
