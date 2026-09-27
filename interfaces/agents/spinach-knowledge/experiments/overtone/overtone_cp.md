# experiments/overtone/overtone_cp.m

- Signature: `spectrum=overtone_cp(spin_system,parameters,H,R,K)`

## Purpose

Cross-polarization overtone experiment. Syntax: spectrum=overtone_cp(spin_system,parameters,H,R,K)

## Physical / mathematical content

The function calculates the overtone reference frequency as `-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`. It projects the supplied `parameters.Nx` and `parameters.Hx` operators over the requested Fokker–Planck spatial dimension.
## Numerical / algorithmic content

The `average` method forms an average pulse Hamiltonian and applies its propagator to `parameters.rho0`; the `fplanck` method applies `shaped_pulse_af` to the initial state. Both methods then delegate frequency-domain acquisition to `overtone_a`.
## Parameters / inputs

- parameters.spins overtone-active nucleus, specified as a single-element cell array
- parameters.spc_dim Fokker-Planck spatial dimension
- parameters.method pulse simulation method, either 'average' or 'fplanck'
- parameters.sweep vector with two elements giving the spectrum frequency extents in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.Nx X Zeeman operator on the quadrupolar nucleus
- parameters.Hx X Zeeman operator on the spin-1/2 nucleus
- parameters.rf_frq spin-lock frequency offset from the overtone frequency on the quadrupolar nucleus, Hz
- parameters.rf_pwr a vector of spin-lock powers on the quadrupolar nucleus (first element) and the spin-1/2 nucleus (second element), rad/s
- parameters.rf_dur spin-lock pulse duration, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the resulting spectrum
- Notes: relaxation must be present in the system dynamics, or the
- matrix inversion in overtone_a function call would fail to
- converge. The relaxation matrix must *not* be thermalised.

## Implementation structure

Validates the inputs, projects the pulse operators, applies the selected spin-lock pulse method, and calls `overtone_a` to acquire the spectrum.