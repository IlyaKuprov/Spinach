# kernel/pulses/grad_sandw.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/grad_sandw.m) · [Spin Dynamics Wiki: grad_sandw.m](https://spindynamics.org/wiki/index.php?title=grad_sandw.m)

Signature: `rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)`

## Purpose and assumptions

Computes the sample-averaged effect of a two-gradient sandwich using Edwards' formalism. It assumes negligible diffusion, linear gradients, and gradients antisymmetric about the sample midpoint. Because the spatial coordinate is integrated out, the returned state is not a spatially resolved state for a later gradient to refocus; this routine is for standalone gradient pairs. Use the imaging context for more sophisticated gradient evolution.

## Inputs and units

- `spin_system`, `L`, and `rho`: a Liouville-space spin system, its Liouvillian, and the state vector. The source accepts the `sphten-liouv` and `zeeman-liouv` formalisms.
- `P`: the source describes this as the total propagator for all events between the two gradients.
- `g_amps`: two real gradient amplitudes in gauss/cm; the source header calls these a row vector, while validation checks for two elements.
- `s_len`: positive sample length in cm.
- `g_durs`: two real, non-negative gradient durations in seconds.
- `s_facs`: two real, non-negative shape factors. Keep the documented `[1 1]` for square gradient pulses.

## What the implementation does

The effective gradient operators use `1e-4 * shape_factor * gradient_amplitude * sample_length * duration * (carrier / magnet)`; the source warns that shifts are ignored. Before propagating, it checks that `L*G_i - G_i*L` has `cheap_norm` no greater than `1e-6` for each effective gradient operator; otherwise it errors. The state is evolved under `L` for each gradient duration, while a block evolution combines the gradient operators with `P` and maps the result back to the state-vector space. The half-gradient evolution factors account for the normalised sample integral.

The routine calls `report` with a progress message. It contains no explicit plotting or file-writing operation.

## References

- [Luke et al., Journal of Magnetic Resonance (2014)](http://dx.doi.org/10.1016/j.jmr.2014.01.011)
- [Spin Dynamics Wiki: `grad_sandw.m`](https://spindynamics.org/wiki/index.php?title=grad_sandw.m)
