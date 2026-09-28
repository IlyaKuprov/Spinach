# kernel/pulses/grad_sandw.m

- Signature: `rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)`

## Purpose

Approximates the effect of two linear gradient pulses with intervening evolution on the sample-averaged state, using Edwards' formalism. It assumes negligible diffusion and gradients antisymmetric about the sample midpoint. The routine integrates over the spatial coordinate and returns a spatially averaged state, so this result cannot be followed by another gradient pulse to refocus the dephasing.

## Inputs

- `spin_system` — Spinach system in Liouville space (`sphten-liouv` or `zeeman-liouv`).
- `L` — numeric system Liouvillian.
- `rho` — numeric state vector or matrix of states.
- `P` — numeric total propagator for events between the two gradients.
- `g_amps` — two real gradient amplitudes in Gauss/cm.
- `s_len` — positive real sample length in cm.
- `g_durs` — two real, non-negative gradient durations in seconds.
- `s_facs` — two real, non-negative gradient shape factors; use `[1 1]` for square pulses.

The effective gradient operators are built from the carrier frequencies, sample length, pulse amplitudes, durations, and shape factors; frequency shifts are ignored. The Liouvillian must commute with each effective gradient operator to within the source's `1e-6` norm check, or the function errors.

## Output

- `rho` — state vector or state matrix integrated over the sample coordinate after the gradient sandwich.

This function is intended for standalone gradient pairs. Use the imaging context for more sophisticated gradient evolution.

## References

- [Luke et al., Journal of Magnetic Resonance (2014)](http://dx.doi.org/10.1016/j.jmr.2014.01.011)
- [Spin Dynamics Wiki: `grad_sandw.m`](https://spindynamics.org/wiki/index.php?title=grad_sandw.m)
