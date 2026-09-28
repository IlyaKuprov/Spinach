# examples/imaging/spin_echo_1d.m

- Signature: `spin_echo_1d()`

## Purpose

Simulates a 1D spin echo under a gradient, with diffusion and flow. The source estimates seconds of runtime.

## Model and sequence

The system is one `1H` at 5.9 T with chemical-shift value 1.0, using the `sphten-liouv` formalism without a basis approximation. The sample is 0.30 m long with 100 points. The gradient amplitude is `5e-6`, each gradient step lasts `2e-4 s`, and there are 200 gradient steps. The initial and detection state phantoms are uniform; the flow field is set to one and diffusion to `1e-6`.

## Computation and output

The echo is calculated with `imaging(spin_system,@spin_echo,parameters)`. The plotted real echo signal uses a time axis from 0 through 400 gradient-step intervals, labelled in seconds; the figure labels the vertical axis as echo intensity in arbitrary units.

## Attribution

Ahmed Allami; Ilya Kuprov.
