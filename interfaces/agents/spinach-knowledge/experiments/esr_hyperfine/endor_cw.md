# experiments/esr_hyperfine/endor_cw.m

- Signature: `fid=endor_cw(spin_system,parameters,H,R,K)`

## Purpose

Fast approximate simulation of isotropic continuous-wave ENDOR. The calculation records an NMR spectrum weighted by hyperfine couplings.

## Parameters / inputs

- `spin_system` — spin system passed to the simulation and operator functions.
- `parameters.sweep` — nuclear frequency sweep width, Hz.
- `parameters.npoints` — number of FID points to be computed.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.

## Outputs

- `fid` — free induction decay whose Fourier transform approximates a CW ENDOR spectrum.

## Implementation structure

The function moves the inputs into the adjoint representation if needed and checks that the formalism is `sphten-liouv` or `zeeman-liouv`, that `H`, `R`, and `K` are matrices of the same dimensions, and that `parameters.sweep` and `parameters.npoints` each have one element. It composes the Liouvillian as `L=H+1i*R+1i*K`.

The initial state is a sum of nuclear `Lz` states weighted by the absolute values of electron–nuclear hyperfine coupling amplitudes, then normalized. A nuclear `Sy` operator applies a `pi/2` pulse. The function detects nuclear `L+` coherence and acquires `parameters.npoints` FID points with evolution time step `1/parameters.sweep`. It returns the real part of the FID for frequency symmetrization.

<https://spindynamics.org/wiki/index.php?title=endor_cw.m>