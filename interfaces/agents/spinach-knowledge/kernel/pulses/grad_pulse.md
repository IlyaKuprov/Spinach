# kernel/pulses/grad_pulse.m

- Signature: `rho=grad_pulse(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)`

## Purpose

Approximates the effect of a single linear gradient pulse on the sample-averaged density matrix using Edwards' formalism. It assumes negligible diffusion and a gradient antisymmetric about the sample midpoint. The function integrates over the sample coordinate, so a later gradient pulse cannot refocus the spatial dephasing produced here.

## Inputs

- `spin_system` — Spinach system structure in Liouville space (`sphten-liouv` or `zeeman-liouv`).
- `L` — numeric system Liouvillian.
- `rho` — numeric state vector or matrix of states.
- `g_amp` — real scalar gradient amplitude in Gauss/cm.
- `s_len` — positive real sample length in cm.
- `g_dur` — non-negative real gradient duration in seconds.
- `s_fac` — non-negative real shape factor; use 1 for a square gradient pulse.

The calculation constructs an effective gradient operator from the carrier frequencies, sample length, amplitude, duration, and shape factor; frequency shifts are ignored. It requires `L` to commute with that operator to within the source's `1e-6` norm check, and errors otherwise.

## Output

- `rho` — state vector or state matrix integrated over the spatial coordinate after the gradient pulse.

For a gradient sandwich or more sophisticated gradient evolution, use `grad_sandw.m` or the imaging context rather than chaining calls to this spatially averaged result.

## References

- [Luke et al., Journal of Magnetic Resonance (2014)](http://dx.doi.org/10.1016/j.jmr.2014.01.011)
- [Spin Dynamics Wiki: `grad_pulse.m`](https://spindynamics.org/wiki/index.php?title=grad_pulse.m)
