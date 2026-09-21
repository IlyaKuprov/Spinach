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
