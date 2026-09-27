# experiments/esr_dipolar/sifter.m

- Signature: `fid=sifter(spin_system,parameters,H,R,K)`

## Purpose

SIFTER pulse sequence. `H` is the Hamiltonian matrix, `R` is the relaxation matrix, and `K` is the chemical kinetics matrix.

## Physical / mathematical content

The sequence applies a 90-degree X pulse, evolves through the first part of `t1`, applies a 180-degree +X pulse, and evolves through the second part of `t1`. It then applies a 90-degree +Y pulse, evolves through the first part of `t2`, applies a 180-degree +X pulse, and detects during the second part of `t2`.

## Numerical / algorithmic content

- The Liouvillian is assembled as `L=H+1i*R+1i*K`.
- Evolution uses `parameters.timestep` and `parameters.npoints/2-1` steps for each period. The first part of `t1` uses `trajectory` mode; the second part of `t1` and first part of `t2` use `refocus` mode. The state-stack columns are reversed before the first part of `t2`.
- The second part of `t2` uses `observable` mode with `parameters.coil` for 2D detection.

## Parameters / inputs

- `spin_system` — passed to the pulse and evolution operations.
- `H` — Hamiltonian matrix.
- `R` — relaxation matrix.
- `K` — chemical kinetics matrix.
- `parameters.npoints` — number of points in time evolution; must be a finite even real integer greater than or equal to 2.
- `parameters.timestep` — simulation time step, seconds; must be a positive real scalar.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.pulse_opx` — pulse operator in X phase.
- `parameters.pulse_opy` — pulse operator in Y phase.

`H`, `R`, and `K` must be numeric matrices of the same dimensions. All listed `parameters` fields are required.

## Outputs

- `fid` — a 2D free induction decay.

## Implementation structure

The function checks input consistency, composes the Liouvillian, then applies the pulses and evolution periods in sequence. The final evolution call returns `fid`.

Source: <https://spindynamics.org/wiki/index.php?title=sifter.m>