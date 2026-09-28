# examples/imaging/gradient_echo_1d.m

- Signature: `gradient_echo_1d()`

## Purpose

Simulates a one-dimensional gradient echo in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami and Ilya Kuprov.

## Physical / mathematical content

The example uses a single `1H` spin at a magnetic induction of 5.9 with a chemical shift of 1.0. It sets `parameters.dims=0.30`, `parameters.npts=100`, and `parameters.deriv={'period',3}`. Diffusion is set to `1e-6`, and the flow profile is one at every spatial point. The initial-state and detection phantoms are uniform across the sample, with `Lz` used for the initial spin state and `L+` for detection. The relaxation phantoms are empty.

## Numerical / algorithmic content

The spin system uses the `sphten-liouv` formalism with no basis approximation. The `grad_echo` sequence is run through `imaging` with a zero offset, gradient amplitude of `5e-6`, step duration of `2e-4`, and 200 gradient steps. The real part of the resulting echo is plotted against time in seconds, from minus to plus the gradient duration.

## Implementation structure

`gradient_echo_1d` creates the spin system and basis, sets the sequence and sample parameters, calls `imaging(spin_system,@grad_echo,parameters)`, and plots the echo intensity.