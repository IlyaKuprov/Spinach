# experiments/esr_dipolar/deer_3p_hard_echo.m

- Signature: `echo=deer_3p_hard_echo(spin_system,parameters,H,R,K)`

## Purpose and sequence

Sample a three-pulse hard-pulse echo window. Starting from `parameters.rho0`, apply the probe `pi/2` pulse and evolve for `parameters.tb`; apply the pump `pi` pulse and evolve for `parameters.ta-parameters.tb`; apply the probe `pi` pulse, evolve for `parameters.ta-parameters.tc/2`, then detect with `parameters.coil` over the window `parameters.tc` using `parameters.nsteps` samples. `H`, `R` and `K` define the supplied evolution; this routine does not construct a dipolar interaction.

## Parameters and output

- `parameters.ex_prob` and `parameters.ex_pump`: probe and pump pulse operators.
- `parameters.ta`, `parameters.tb` and `parameters.tc`: sequence timing values in seconds.
- `parameters.rho0`: initial state; `parameters.nsteps`: number of points across the detection window; `parameters.coil`: detection operator.
- Returns `echo`, the sampled observable across the echo window.

## Input requirements

`H`, `R` and `K` must be same-sized matrices. `parameters.ta`, `parameters.tb` and `parameters.tc` must be positive, `parameters.tb <= parameters.ta`, `parameters.tc/2 <= parameters.ta`, and `parameters.nsteps` must be a positive integer.
