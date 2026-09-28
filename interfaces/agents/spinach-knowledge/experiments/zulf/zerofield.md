# experiments/zulf/zerofield.m

- Signature: `fid=zerofield(spin_system,parameters,H,R,K)`

## Purpose

Simulates a Budker-group-style gamma-weighted pulse-acquire sequence in zero field. The initial state, pulse operator, and detection state are weighted by each spin's magnetogyric ratio relative to 1H, as in preparation with a high-field pre-polarisation magnet.

## Physical / mathematical content

The supplied matrices form the Liouvillian `L = H + 1i*R + 1i*K`. The initial state is the weighted sum of single-spin `Lz` states, while the coil state is the weighted sum of `L+` states. A weighted `Ly` pulse operator acts on the initial state with the requested flip angle. Detection is performed against the weighted coil state.

## Numerical / algorithmic content

The acquisition time step is `1/parameters.sweep`; the evolution call uses `parameters.npoints-1` propagation steps and observable mode. The routine requires positive scalar sweep width, a positive integer point count, a real scalar flip angle, and a detection mode of `'uniaxial'` or `'quadrature'`. The supplied `H`, `R`, and `K` must be numeric matrices of equal size. In quadrature mode the complex signal is retained; uniaxial detection takes its real part.

## Parameters / inputs

- `parameters.sweep` - the width of the spectral window (Hz)
- `parameters.npoints` - number time steps in the simulation
- `parameters.detection` - 'uniaxial' to emulate common ZULF hardware, 'quadrature' for proper frequency sign discrimination
- `parameters.flip_angle` - pulse flip angle in radians for protons; for other nuclei, this will be scaled by the gamma ratio
- `H` - Hamiltonian matrix, received from context function
- `R` - relaxation superoperator, received from context function
- `K` - kinetics superoperator, received from context function

## Outputs

- `fid` - free induction decay

## Implementation structure

After consistency checks, the routine composes the Liouvillian, constructs the gamma-weighted initial state and detection coil, builds the weighted `Ly` pulse operator, applies the pulse, and calls `evolution` for acquisition. The detection-mode switch leaves quadrature data unchanged and removes the imaginary part for uniaxial detection.
