# experiments/esr_dipolar/oopeseem.m

- MATLAB implementation: [experiments/esr_dipolar/oopeseem.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/oopeseem.m)

Source: https://spindynamics.org/wiki/index.php?title=oopeseem.m

`fid=oopeseem(spin_system,parameters,H,R,K)`

## What it calculates

This routine simulates an ideal-pulse OOP-ESEEM echo and returns a time-domain signal. ESEEM-family modulation can carry nuclear-frequency information arising from electron–nuclear hyperfine coupling. ENDOR is a different experiment and is not implemented here; the source does not expand the acronym OOP or assign a numerical hyperfine coupling.

The code evaluates the supplied spin system and state. It does not implement field sweeping, DNP/hyperpolarisation, imaging, or measurement acquisition.

## Inputs

- `spin_system` — Spinach spin system; only `sphten-liouv` and `zeeman-liouv` formalisms are accepted.
- `parameters.npoints` — number of computed points; the source checks for one element.
- `parameters.timestep` — simulation time step in seconds; the source checks for one element.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.screen` — optional screen state, documented as the Hermitian conjugate of the detection state; defaults to `[]`.
- `parameters.pulse_op` — caller-supplied pulse operator, dimension-matched to `H`.
- `H`, `R`, `K` — dimension-matched numeric Hamiltonian, relaxation, and kinetics matrices. The function calls `sim2liouv` as needed, then forms `L=H+1i*R+1i*K`.

## Sequence and propagation

The supplied pulse operator is applied to `rho0` with a `pi/4` rotation. The state evolves with `timestep/2` for `npoints-1` steps in trajectory mode with `screen`; a `pi` rotation follows, then a second `timestep/2` evolution for `npoints-1` steps in refocus mode with `coil` as observable. Detection returns the transposed full projection `coil'*rho_stack` as `fid`. The source header describes the ideal pulse propagators as `exp(-i*A*pi/4)` and `exp(-i*A*pi)`.

## Output and scope

`fid` is the time-domain OOP-ESEEM trace; the source does not return a separate time vector or name array axes. Its source-backed pulse-angle values are `pi/4` and `pi` radians; it supplies no concrete experimental parameter values or DOI.

The source header shows `fid=eseem(...)` as its syntax example, whereas the declared function is `fid=oopeseem(...)`.
