# experiments/cp_contact_soft.m

- Signature: `contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)`

## Purpose

Simulates a rotating-frame cross-polarisation contact curve with a soft `pi/2` high-gamma excitation pulse, followed by spin-lock evolution.

## Implementation

The routine composes `L=H+1i*R+1i*K`, wipes the low-gamma spin state from `parameters.rho0`, applies the high-gamma excitation pulse, and evolves during the CP contact while detecting on `parameters.coil`. The requested contact evolution is set by the time step and number of steps.

## Parameters / inputs

- `parameters.spins`: working spins in a cell array, high-gamma first and low-gamma last (for example, {'1H','13C'}).
- `parameters.hi_pwr`: high-gamma excitation-pulse nutation frequency, Hz.
- `parameters.cp_pwr`: two-channel nutation frequencies during CP contact, Hz.
- `parameters.timestep`: CP contact time step, s.
- `parameters.nsteps`: number of CP contact time steps.
- `parameters.rho0`: initial state; the low-gamma spin state is wiped before the sequence.
- `parameters.coil`: detection state vector.
- `H`: Hamiltonian matrix supplied by the context function.
- `R`: relaxation superoperator supplied by the context function.
- `K`: kinetics superoperator supplied by the context function.

## Output

- `contact_curve`: signal detected on the coil state during the CP contact.

[Source page](https://spindynamics.org/wiki/index.php?title=cp_contact_soft.m)
