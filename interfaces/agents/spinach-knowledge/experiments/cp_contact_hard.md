# experiments/cp_contact_hard.m

- Signature: `contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)`

## Purpose

Computes a rotating-frame cross-polarisation contact curve using an ideal pi/2 excitation pulse and hard spin-lock terms.

## Implementation

The routine composes `L=H+1i*R+1i*K`, applies the specified excitation operator(s) to `parameters.rho0` with a pi/2 flip angle, and records the initial coil signal. At each entry of `parameters.time_steps`, it forms the spin-lock contribution from the channel operators and nutation frequencies, propagates the state for that interval, and records the coil-detected signal.

## Parameters / inputs

- `parameters.irr_powers`: matrix of spin-lock nutation frequencies, with channels in rows and time slices in columns, Hz.
- `parameters.irr_opers`: cell array of the spin operators for the spin lock on each channel.
- `parameters.exc_opers`: cell array of spin operators for the ideal pi/2 excitation pulse; the same flip angle is used on all channels.
- `parameters.time_steps`: vector of time-slice durations, s.
- `parameters.rho0`: initial state vector.
- `parameters.coil`: detection state vector.
- `H`: Hamiltonian matrix supplied by the context function.
- `R`: relaxation superoperator supplied by the context function.
- `K`: kinetics superoperator supplied by the context function.

## Output

- `contact_curve`: coil-detected contact curve, including the initial signal and one sample after each time slice.
