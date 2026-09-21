# experiments/zulf/zerofield.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/zulf/zerofield.m`
- Signature: `fid=zerofield(spin_system,parameters,H,R,K)`
- Total lines: 133

## Purpose

Budker group style gamma-weighted pulse-acquire sequence in zero field. Uses gamma-weighted initial state (corresponding to using a pre-polarisation magnet at high temperature), gamma-weighted pulse operators, and gamma-weighted detection state. Syntax: fid=zerofield(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Zero-field experiment implementations. They propagate J-coupled systems in the absence of strong carrier terms and often model abrupt field switching.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.sweep -the width of the spectral window (Hz)
- parameters.npoints -number time steps in the simulation
- parameters.detection -'uniaxial' to emulate common ZULF
- hardware, 'quadrature' for proper
- frequency sign discrimination
- parameters.flip_angle -pulse flip angle in radians for
- protons; for other nuclei, this
- will be scaled by the gamma ratio
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay

## Implementation structure

- Budker group style gamma-weighted pulse-acquire sequence in zero
- field. Uses gamma-weighted initial state (corresponding to using
- a pre-polarisation magnet at high temperature), gamma-weighted
- pulse operators, and gamma-weighted detection state. Syntax:
- fid=zerofield(spin_system,parameters,H,R,K)
- parameters.sweep -the width of the spectral window (Hz)
- parameters.npoints -number time steps in the simulation
- parameters.detection -'uniaxial' to emulate common ZULF
- hardware, 'quadrature' for proper
- frequency sign discrimination
- parameters.flip_angle -pulse flip angle in radians for
- protons; for other nuclei, this

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `weights()`, `state()`, `operator()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `ismember()`, `elseif()`.
