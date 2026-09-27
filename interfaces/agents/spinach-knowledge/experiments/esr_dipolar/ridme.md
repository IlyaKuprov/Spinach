# experiments/esr_dipolar/ridme.m

- Signature: `answer=ridme(spin_system,parameters,H,R,K)`

## Purpose

RIDME pulse sequence using idealized hard pulses that affect only the user-specified electron. `H` is the Hamiltonian commutation superoperator, `R` is the relaxation superoperator, and `K` is the chemical kinetics superoperator.

## Physical / mathematical content

The sequence applies five pulses to the probe spin, with evolution during `tau 1`, `tau 2`, and the mixing time. Phase cycling at the third and fourth pulses produces four signal components. Relaxation must be present for this experiment to work.

## Numerical / algorithmic content

The Liouvillian is `L=H+1i*R+1i*K`. Evolution uses the specified step size and step counts to generate a trajectory and refocus the four phase-cycle branches. The real and imaginary signal components are projections onto the corresponding probe-spin coil states, normalized by their norms.

## Parameters / inputs

- `H`: Hamiltonian, received from the context function.
- `R`: relaxation superoperator, received from the context function.
- `K`: kinetics superoperator, received from the context function.
- `parameters.rho0`: initial state.
- `parameters.probe_spin`: number of the spin on which the sequence operates; a positive real integer.
- `parameters.stepsize`: step size for the increment of the relaxation period, seconds; a positive real scalar.
- `parameters.nsteps(1)`: number of steps for `tau 1`.
- `parameters.nsteps(2)`: number of steps for `tau 2`. `parameters.nsteps` must be a row vector of two positive integers.
- `parameters.tmix`: mixing time, seconds; a positive real scalar.

`H`, `R`, and `K` must be numeric matrices of the same dimensions.

## Outputs

Each of the following has `.real` and `.imag` quadrature components of the signal, corresponding to phase-cycle instances on the third, fourth, and fifth pulses in the RIDME sequence:

- `answer.pxpxpx`
- `answer.pypypx`
- `answer.mxmxpx`
- `answer.mymypx`

## Implementation structure

A `+pi/2` pulse about `Sx` is followed by evolution for `tau 1`, a `+pi` pulse about `Sx`, and trajectory evolution over `parameters.nsteps(1)+parameters.nsteps(2)` steps. The third pulse creates `+Sx`, `+Sy`, `-Sx`, and `-Sy` phase-cycle branches, each of which evolves for `parameters.tmix`. The fourth pulse applies the corresponding `+Sx`, `+Sy`, `-Sx`, or `-Sy` rotation. After refocusing evolution, the fifth pulse applies `+pi` about `Sx` to every branch, followed by evolution for `tau 2` and observation.

<https://spindynamics.org/wiki/index.php?title=ridme.m>