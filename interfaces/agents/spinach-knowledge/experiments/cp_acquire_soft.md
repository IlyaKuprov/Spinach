# experiments/cp_acquire_soft.m

- Signature: `fid=cp_acquire_soft(spin_system,parameters,H,R,K)`

## Purpose

Simulates rotating-frame cross-polarisation followed by time-domain FID acquisition. The source describes wiping the low-gamma spin state before the CP stage and decoupling the high-gamma spins during acquisition.

## Implementation

The routine composes `L=H+1i*R+1i*K`, wipes the low-gamma part of `parameters.rho0`, applies the high-gamma excitation pulse, and evolves during the CP contact. It then acquires the FID on `parameters.coil` with the specified sweep width and point count.

## Parameters / inputs

- `parameters.spins`: working spins in a cell array, high-gamma first and low-gamma last (for example, {'1H','13C'}).
- `parameters.hi_pwr`: high-gamma excitation-pulse nutation frequency, Hz.
- `parameters.cp_pwr`: two-channel nutation frequencies during CP contact, Hz.
- `parameters.cp_dur`: CP contact duration, s.
- `parameters.rho0`: initial state; the low-gamma spin state is wiped before the sequence.
- `parameters.coil`: detection state.
- `parameters.sweep`: FID sweep width, Hz.
- `parameters.npoints`: number of FID points.
- `H`: Hamiltonian matrix supplied by the context function.
- `R`: relaxation superoperator supplied by the context function.
- `K`: kinetics superoperator supplied by the context function.

## Output

- `fid`: signal detected on the coil state during the sequence.
