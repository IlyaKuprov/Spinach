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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 41-42: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 44-45: Magnetogyric ratio weights relative to 1H; implemented by `weights=spin_system.inter.gammas/spin('1H')`.
- Lines 47-48: Get gamma-weighted initial state; implemented by `rho=sparse(0)`.
- Lines 53-54: Get gamma-weighted detection state; implemented by `coil=sparse(0)`.
- Lines 59-60: Get gamma-weighted pulse operator; implemented by `Sy=sparse(0)`.
- Lines 65-66: Apply the pulse; implemented by `rho=step(spin_system,Sy,rho,parameters.flip_angle)`.
- Lines 68-69: Compute the digitisation parameters; implemented by `timestep=1/parameters.sweep`.
- Lines 71-73: Run the simulation; implemented by `fid=evolution(spin_system,L,coil,rho,timestep, parameters.npoints-1,'observable')`.
- Lines 75-76: Emulate detection hardware; implemented by `switch parameters.detection`.
- Lines 80-82: Do nothing; implemented by `case 'uniaxial'`.
- Lines 84-85: Destroy imaginary part; implemented by `fid=real(fid)`.

### Control flow inferred from the code

- Line 49: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 55: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 61: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 76: dispatches on `parameters.detection`; cases `'quadrature'`, `'uniaxial'`.

### Key state/data transformations

- Lines 42: computes `L` using `L=H+1i*R+1i*K`.
- Lines 45: computes `weights` using `weights=spin_system.inter.gammas/spin('1H')`.
- Lines 48: computes `rho` using `rho=sparse(0)`.
- Lines 54: computes `coil` using `coil=sparse(0)`.
- Lines 60: computes `Sy` using `Sy=sparse(0)`.
- Lines 69: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 72-73: computes `fid` using `fid=evolution(spin_system,L,coil,rho,timestep, parameters.npoints-1,'observable')`.

### Local helper functions

- Line 92: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

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
