# experiments/esr_hyperfine/hyscore.m

- Signature: `fid=hyscore(spin_system,parameters,H,R,K)`

## Purpose

HYSCORE experiment, implemented as described in Szosenfogel and Goldfarb (http://dx.doi.org/10.1080/00268979809483260). Syntax: fid=hyscore(spin_system,parameters,H,R,K)

## Physical / mathematical content


- Implements a HYSCORE pulse sequence with electron pi/2 rotations about x, a tau delay, zero-electron-coherence filtering, indirect evolution, and a final electron pi pulse.
- The evolution generator is `L = H + iR + iK`; the detection state uses the coil operator propagated backward through the detection interval and a negative pi/2 pulse.

## Numerical / algorithmic content


- The two indirect dimensions are sampled on the configured time grids; the function returns a two-dimensional time-domain FID.
- Uses `step()`, `coherence()`, and `evolution()` for pulse and coherence selection; Fourier transformation and spectrum post-processing are not performed here.

## Parameters / inputs

- parameters.nsteps number of points to be computed
- in each dimension
- parameters.sweep sweep width, Hz
- parameters.tau tau delay, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay that Fourier
- transforms into a HYSCORE spectrum
- Note: the sequence uses ideal pulses, replace with shaped_pulse_af()
- to have soft pulses instead.

## Implementation structure


- Converts the input to the adjoint representation when needed, builds the electron rotation and coherence-selection operators, and validates the Liouville-space dimensions and `parameters.nsteps`.
- Executes the two indirect evolution periods and evaluates the detection trajectory to populate the HYSCORE FID.
