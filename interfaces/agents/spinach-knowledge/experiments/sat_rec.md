# experiments/sat_rec.m

- Signature: `fids=sat_rec(spin_system,parameters,H,R,K)`

## Purpose

Computes a saturation-recovery pulse sequence with analytical saturation, using the unit state as the initial condition.

## Numerical / algorithmic content

The routine composes `L=H+1i*R+1i*K`, starts from `unit_state(spin_system)`, and propagates a trajectory over `parameters.n_delays` equally spaced relaxation periods spanning `parameters.max_delay`. It applies a 90-degree pulse about the source-defined `Ly` operator to each trajectory state, then acquires an FID for each state using the detection state on `parameters.spins{1}`, dwell time `1/parameters.sweep`, and `parameters.npoints-1` intervals. The FIDs are returned as columns.

## Parameters / inputs

- `parameters.sweep` — spectrum sweep width, Hz
- `parameters.npoints` — number of points in each FID
- `parameters.spins` — nuclei on which the sequence runs, specified as {'1H'}, {'13C'}, etc.
- `parameters.max_delay` — longest relaxation delay
- `parameters.n_delays` — number of relaxation delays to run
- `H` — Hamiltonian matrix, received from the context function
- `R` — relaxation superoperator, received from the context function; it must be thermalised
- `K` — kinetics superoperator, received from the context function

## Outputs

- `fids` — free induction decays for each delay starting from zero, with individual FIDs in columns

## Credit

Zak El-Machachi

## Reference

- [Spinach documentation](https://spindynamics.org/wiki/index.php?title=sat_rec.m)
