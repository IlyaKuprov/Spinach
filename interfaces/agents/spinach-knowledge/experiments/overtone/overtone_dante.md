# experiments/overtone/overtone_dante.m

- Signature: `spectrum=overtone_dante(spin_system,parameters,H,R,K)`

## Purpose

Overtone DANTE experiment with frequency-domain acquisition.

## Physical / mathematical content

The routine calculates the overtone reference frequency as `-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`, forms `L=H+1i*R+1i*K`, and extends `parameters.Lx` across `parameters.spc_dim` spatial dimensions.
## Numerical / algorithmic content

The function checks its inputs and rejects a pulse that does not fit in a rotor cycle. It constructs pulse and free-evolution propagators, combines them into one cycle, applies that cycle `parameters.n_periods*parameters.pulse_num` times with `multiprop`, and calls `overtone_a` for acquisition.
## Syntax

```matlab
spectrum=overtone_dante(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.pulse_dur -duration of the pulse, seconds
- parameters.pulse_amp -amplitude of the pulse, rad/s
- parameters.pulse_num -number of pulses within the rotor period
- parameters.n_periods -number of rotor periods that the sequence is active for
- parameters.spins -overtone-active nucleus, specified as a single-element cell array
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.Lx -X Zeeman operator on the quadrupolar nucleus
- parameters.rf_frq -pulse frequency offset from the overtone frequency, Hz
- parameters.rate -rotor frequency in Hz
- parameters.sweep -acquisition sweep range, Hz
- parameters.npoints -number of acquisition points
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -overtone spectrum

## Implementation structure

After validation, the routine builds the Liouvillian and spatially extended pulse operator, constructs the DANTE pulse train from the pulse and free-evolution propagators, propagates `parameters.rho0`, and delegates spectral acquisition to `overtone_a`.