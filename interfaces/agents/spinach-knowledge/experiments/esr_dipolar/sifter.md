# experiments/esr_dipolar/sifter.m

- MATLAB implementation: [experiments/esr_dipolar/sifter.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/sifter.m)

Source: https://spindynamics.org/wiki/index.php?title=sifter.m

`fid=sifter(spin_system,parameters,H,R,K)`

## What it calculates

This function simulates a SIFTER pulse sequence and returns a two-dimensional free-induction decay (FID). It uses the supplied spin-system Hamiltonian and pulse operators; the source does not define a field sweep, DNP/hyperpolarisation step, imaging dimension, or a measured acquisition. No acronym expansion or specific hyperfine coupling is given in the routine source.

## Inputs

- `spin_system` — Spinach spin system.
- `parameters.npoints` — number of points; it must be a finite even real integer of at least 2.
- `parameters.timestep` — positive real time step in seconds.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.pulse_opx` and `parameters.pulse_opy` — caller-supplied X- and Y-phase pulse operators, dimension-matched to `H`.
- `H`, `R`, `K` — dimension-matched numeric matrices for the Hamiltonian, relaxation, and chemical kinetics contributions. The generator is `L=H+1i*R+1i*K`.

## Sequence and propagation

The routine applies a `pi/2` X pulse, evolves a first interval with `npoints/2-1` steps at `timestep` in trajectory mode, and applies a `pi` X pulse. It then refocus-evolves the stored stack for the first part of the echo, applies a `pi/2` Y pulse, reverses the stored order, and refocus-evolves the second delay period. A final `pi` X pulse precedes observable-mode evolution with `coil` over `npoints/2-1` steps. Pulse angles are in radians; the user supplies the X and Y pulse operators.

## Output and scope

`fid` is the source-documented 2D FID. Its two sampling dimensions arise from the two evolution periods, but the function does not return separate named axis vectors; `timestep` and the even `npoints` control sampling. The numeric examples specified by the routine are the pulse angles `pi/2` and `pi`, and `npoints` must be at least 2 and even. No DOI or concrete physical parameter set is supplied.
