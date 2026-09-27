# experiments/slowpass.m

- Signature: `spectrum=slowpass(spin_system,parameters,H,R,K)`

## Purpose

Calculates spectrum values at the frequency positions specified by `parameters.sweep`, without first calculating the complete free induction decay.

## Numerical / algorithmic content

The routine constructs a frequency grid, projects the initial state, detection state, and Liouvillian into each selected subspace, and solves a linear system at each frequency. Depending on the configured execution path it uses backslash or preconditioned GMRES; a GPU backslash path is also provided. The accumulated spectrum is scaled by the sampling rate implied by the frequency grid to match the unnormalised FFT amplitude convention. The relaxation matrix `R` must not be thermalised.

## Parameters / inputs

- `parameters.sweep` — two-element vector giving the spectrum frequency extents, Hz
- `parameters.npoints` — number of points in the spectrum
- `parameters.rho0` — initial state
- `parameters.coil` — detection state
- `H` — Hamiltonian matrix, received from the context function
- `R` — relaxation superoperator, received from the context function; it must not be thermalised
- `K` — kinetics superoperator, received from the context function

## Outputs

- `spectrum` — spectrum of the system for the specified starting state and detection state over the requested frequency interval

## Reference

- [Spinach documentation](https://spindynamics.org/wiki/index.php?title=slowpass.m)
