# examples/optimal_control/features_dt_var.m

- Signature: `features_dt_var()`
- Source: [examples/optimal_control/features_dt_var.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_dt_var.m)

## Purpose and model

This simulated optimal-control example transfers normalised longitudinal magnetisation from 1H to 19F through the scalar-coupled 1H-13C-19F model. It sets a 9.4 T field, zero chemical shifts, and couplings of 140 Hz (1H-13C) and -160 Hz (13C-19F). The model is defined in the script; it does not import measured data.

## Controls and optimisation

The six controls are x and y RF components on each nucleus. The pulse has 50 nonuniform intervals, with durations defined by `3e-4*(0.25+0.75*cos(linspace(-pi/2,pi/2,50)))` seconds (about 75 to 300 microseconds). Optimisation uses five RF-power levels from `2*pi*800` to `2*pi*1200` rad/s, the SNS penalty with weight 100, the `lbfgs` method, and a 100-iteration limit. A random 6-by-50 initial waveform is passed to `fmaxnewton` with `@grape_xy`; correlation-order, per-spin, and x/y-control plots are enabled.

## Output and limits

The optimised waveform is scaled by the mean power level and propagated once with `shaped_pulse_xy` using `expv-pwc`. The reported quantity is the real overlap `Re[rho_targ'*rho(T)]`. This final check uses a single drift Hamiltonian and the mean-power-scaled pulse; it is not a reported sweep over all five powers. The example is a model simulation, not a hardware or experimental validation.
