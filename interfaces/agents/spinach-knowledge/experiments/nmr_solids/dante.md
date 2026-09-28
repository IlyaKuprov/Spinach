# experiments/nmr_solids/dante.m

- Signature: `fid=dante(spin_system,parameters,H,R,K)`

## Purpose

`dante` implements a DANTE pulse sequence and is normally called through the `singlerot.m` context, which supplies `H`, `R`, and `K`.

## Physical / mathematical content

The sequence acts on the configured working spin, with decoupling isotopes supplied separately. Pulse timing and repetition are set by the pulse parameters, rotor frequency, and number of rotor periods.

## Numerical / algorithmic content

The source applies timed pulses with free evolution between them and records the acquisition signal. It does not implement powder-orientation averaging or cross-polarisation/recoupling.

## Parameters / inputs

- parameters.pulse_dur -duration of each pulse, seconds
- parameters.pulse_amp -amplitude of each pulse, Hz
- parameters.pulse_num -number of pulses within rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -working spin, specified as a
- single-element cell array
- parameters.decouple -isotopes to decouple, specified
- as a cell array
- parameters.rate -rotor frequency in Hz
- parameters.sweep -acquisition sweep width in Hz
- parameters.npoints -number of acquisition points
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+

## Outputs

- fid -free induction decay
