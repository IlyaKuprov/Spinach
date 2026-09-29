# experiments/esr_dipolar/ridme.m

- MATLAB implementation: [experiments/esr_dipolar/ridme.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/ridme.m)

Source: https://spindynamics.org/wiki/index.php?title=ridme.m

`answer=ridme(spin_system,parameters,H,R,K)`

## What it calculates

RIDME is a relaxation-induced dipolar-modulation experiment. This routine implements an ideal hard-pulse, five-pulse phase-cycled sequence on the selected probe spin. Its source explicitly notes that relaxation must be present for the experiment to work. This is not an ENDOR or DNP/hyperpolarisation workflow, a field sweep, or an imaging routine.

## Inputs

- `spin_system` — Spinach spin system.
- `parameters.rho0` — initial state.
- `parameters.probe_spin` — positive integer spin index for the observed spin.
- `parameters.stepsize` — positive time increment in seconds.
- `parameters.nsteps` — row vector of two positive integers; entries count steps for the `tau1` and `tau2` segments of one combined delay sweep.
- `parameters.tmix` — positive mixing time in seconds.
- `H`, `R`, `K` — dimension-matched numeric generator matrices for the Hamiltonian, relaxation, and chemical kinetics contributions; the routine forms `L=H+1i*R+1i*K`. Relaxation is part of the required physical model, not an optional post-processing correction.

## Sequence and propagation

The implementation constructs the probe-spin X/Y pulse operators and real/imaginary detection quadratures. It applies a `pi/2` X pulse, evolves for `stepsize*nsteps(1)`, applies a `pi` X pulse, and samples one trajectory with `nsteps(1)+nsteps(2)` steps of duration `stepsize`. Four third-pulse branches use `+X`, `+Y`, `-X`, and `-Y` rotations of `pi/2`; each evolves for `tmix` and receives its corresponding fourth pulse. The stored evolution is reversed for refocusing, followed by a `pi` X pulse and final evolution over `stepsize*nsteps(2)`. These pulse rotations are in radians.

## Output and limitations

The result has four phase-cycle channels: `answer.pxpxpx`, `answer.pypypx`, `answer.mxmxpx`, and `answer.mymypx`; each has `real` and `imag` quadrature components. The code projects and normalises by the corresponding probe-spin coil-state norm. Every real/imaginary component is a row vector of length `nsteps(1)+nsteps(2)+1`, including the initial trajectory point; no separate named delay-axis vectors are returned. The source supplies pulse angles and parameter constraints, but no example numerical times, measured data, or DOI.
