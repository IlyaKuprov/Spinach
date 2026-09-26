# kernel/pulses/isergen.m

- Signature: `H=isergen(HL,HM,HR,dt)`

## Purpose

2nd and 4th order Iserles product quadrature generators for one time propagation step in the case of state-inde- pendent Hamiltonian. Syntax: H=isergen(HL,HM,HR,dt)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
