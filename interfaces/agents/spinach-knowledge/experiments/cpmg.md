# experiments/cpmg.m

- Signature: `fid=cpmg(spin_system,parameters,H,R,K)`

## Purpose

Simulates a CPMG echo train and detects the signal throughout the sequence.

## Implementation

The routine composes `L=H+1i*R+1i*K`, propagates an initial half-echo, and records the coil signal. It then repeats the requested number of CPMG loops: apply `parameters.pulse_op` with a pi rotation, propagate the next half-echo, and append the detected signal. The half-echo duration is represented by `parameters.timestep` and `parameters.npoints`.

## Parameters / inputs

- `parameters.rho0`: initial state.
- `parameters.coil`: detection state.
- `parameters.pulse_op`: pulse operator used for the CPMG refocusing pulses.
- `parameters.nloops`: number of CPMG loops.
- `parameters.timestep`: time step.
- `parameters.npoints`: number of propagation steps per half-echo.
- `H`: Hamiltonian matrix supplied by the context function.
- `R`: relaxation superoperator supplied by the context function.
- `K`: kinetics superoperator supplied by the context function.

## Output

- `fid`: free induction decay recorded throughout the sequence.

[Source page](https://spindynamics.org/wiki/index.php?title=cpmg.m)
