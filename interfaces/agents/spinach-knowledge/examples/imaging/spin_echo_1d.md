# examples/imaging/spin_echo_1d.m

- Signature: `spin_echo_1d()`
- Source: [examples/imaging/spin_echo_1d.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/spin_echo_1d.m)

## Experiment

This driver simulates a one-dimensional spin echo under a gradient, with diffusion and flow. The source comment estimates a runtime of seconds.

## Model and sequence parameters

The system is a single 1H with magnetic-induction parameter 5.9 and chemical-shift value 1.0; the source does not label units for these values. The basis uses sphten-liouv with no approximation. The source does not specify a T1/T2 relaxation model for this example.

The spatial model has dims=0.30 and 100 sample points; dims has no unit annotation in the source. The gradient amplitude is 5e-6, applied for 200 steps of 2e-4 seconds each. The gradient-amplitude unit is not stated. The spatial derivative option is {'period',3}. Initial magnetisation and receiver sensitivity use uniform phantoms with Lz preparation and L+ detection operators. The relaxation phantom is zero, with relaxation(spin_system) supplied as its operator. Flow is set to 1 at every point and diffusion to 1e-6; the driver does not state their units.

## Computation and observable

The driver invokes imaging(spin_system,@spin_echo,parameters) and plots the real echo signal. Its plot constructs a 401-point time axis, `(0:400)*g_step_dur`, which spans 0 to 0.08 seconds, and labels the vertical axis as echo intensity in arbitrary units. This describes the requested plot, not a reported numerical echo trace; no simulation was run for this page.

## Scope

The example is a one-dimensional spin-echo acquisition with gradient, flow, and diffusion settings. It does not state a measured benchmark or a comparison with a reference echo.

## Attribution

Ahmed Allami; Ilya Kuprov.
