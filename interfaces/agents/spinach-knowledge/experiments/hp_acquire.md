# experiments/hp_acquire.m

- Signature: `fid=hp_acquire(spin_system,parameters,H,R,K)`

## Purpose

Standard pulse-acquire sequence with a hard pulse. The user supplies the pulse operator, pulse angle and initial condition. Echo detection is optional.

## Physical / mathematical content

The Liouvillian is assembled as `L=H+1i*R+1i*K`. The hard pulse acts on the initial state; when echo detection is requested, the sequence is `echo_time - pulse - echo_time - fid`.

## Numerical / algorithmic content

The function checks input consistency, projects the pulse operator using `kron(speye(parameters.spc_dim),parameters.pulse_op)`, and applies the pulse with `step`. If `parameters.echo_time` is supplied, it evolves for the echo time, projects and applies the echo pulse, then evolves for the echo time again. It then applies decoupling and records the coil observable using `evolution` with a sampling interval of `1/parameters.sweep` and `parameters.npoints-1` evolution steps.

## Parameters / inputs

- `parameters.sweep` — sweep width, Hz; a positive real scalar.
- `parameters.npoints` — number of points in the FID; a positive integer.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.pulse_op` — pulse operator.
- `parameters.pulse_angle` — pulse angle in radians; a real scalar.
- `parameters.decouple` — spins to decouple, e.g. `{'15N','13C'}`; a cell array of isotope strings, or an empty cell array. Nonempty analytical decoupling is available only in the `sphten-liouv` formalism.
- `parameters.echo_time` — optional positive real echo time for echo detection (`echo_time - pulse - echo_time - fid`).
- `parameters.echo_oper` — optional pulse operator for echo detection; required when `parameters.echo_time` is supplied.
- `parameters.echo_angle` — optional real pulse angle for echo detection; required when `parameters.echo_time` is supplied.
- `H` — Hamiltonian matrix received from the context function.
- `R` — relaxation superoperator received from the context function.
- `K` — kinetics superoperator received from the context function.

`H`, `R` and `K` must be numeric matrices of the same size. If either echo pulse parameter is supplied, `parameters.echo_time` must also be supplied.

## Outputs

- `fid` — free induction decay observed through the detection state specified in `parameters.coil`.

## Implementation structure

Input validation precedes Liouvillian construction, pulse application, optional echo evolution, decoupling and observable acquisition.

## Authors and link

- ilya.kuprov@weizmann.ac.il
- ledwards@cbs.mpg.de
- <https://spindynamics.org/wiki/index.php?title=hp_acquire.m>