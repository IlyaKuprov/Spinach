# experiments/esr_dipolar/deer_3p_hard_deer.m

- Signature: `deer=deer_3p_hard_deer(spin_system,parameters,H,R,K)`

## Purpose and sequence

Generate a three-pulse DEER trace using the spin system, state, pulse operators and evolution matrices supplied by the caller. Starting from `parameters.rho0`, apply a probe `pi/2` pulse, evolve, apply a pump `pi` pulse through the pump-pulse sandwich, evolve, apply the probe `pi` pulse, and detect on the probe spin with `parameters.coil_prob`. Evolution uses the supplied `H`, `R` and `K` matrices; the function does not construct a dipolar interaction. Hard pulses are appropriate only for spin-1/2 systems; for higher-spin systems, supply transition-selective pulse operators.

## Parameters and output

- `parameters.ex_prob` and `parameters.ex_pump`: probe and pump excitation operators.
- `parameters.stepsize` and `parameters.nsteps`: time increment and number of steps in the pump-pulse sandwich; the increment must be positive and the count a positive integer.
- `parameters.output`: choose `brief` or `detailed`.
- `brief` returns `deer.deer_trace`. `detailed` additionally requires `parameters.ex_hard`, `parameters.coil_pump`, `parameters.spectrum_sweep` and `parameters.spectrum_nsteps`, and returns `deer.hard_pulse_fid`, `deer.prob_pulse_fid` and `deer.pump_pulse_fid` as well.

## Input requirements

`H`, `R` and `K` must be same-sized matrices. Supply the initial state and probe detection state as `parameters.rho0` and `parameters.coil_prob`.
