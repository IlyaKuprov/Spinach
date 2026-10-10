# kernel/pulses/grad_pulse.m

[Source: `kernel/pulses/grad_pulse.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/grad_pulse.m)

- Signature: `rho=grad_pulse(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)`

## Purpose

Approximates the effect of one linear gradient pulse on the sample-averaged density matrix using Edwards' formalism. It assumes negligible diffusion and a gradient antisymmetric about the sample midpoint. Because the returned state has already been averaged over sample position, it cannot retain spatial information needed to model subsequent gradient refocusing; use `grad_sandw.m` or the imaging context for that work.

## Inputs and timing

- `spin_system` — Spinach system structure in Liouville space.
- `L` — system Liouvillian; `rho` — state vector.
- `g_amp` — gradient amplitude in Gauss/cm; `s_len` — sample length in cm.
- `g_dur` — gradient duration in seconds; `s_fac` — non-negative gradient shape factor, with 1 for a square gradient.

The implementation first propagates the state under `L` for `g_dur`. It forms a gradient operator proportional to `g_amp*s_len*g_dur*s_fac`, then integrates the spatial coordinate over the normalised sample interval using an auxiliary block-operator evolution. This is an analytic sample average in the Edwards approximation, not a discretised user-supplied gradient waveform. The source notes that chemical shifts are ignored in the gradient operator.

The function accepts numeric arguments apart from `spin_system`, requires a Liouville-space formalism, and checks the gradient amplitude is real scalar, sample length positive, and duration and shape factor non-negative scalars. It also requires the Liouvillian and gradient operator to commute within the source's `1e-6` norm threshold; otherwise it errors. Progress is reported through Spinach's `report` routine; no file is written.

## References

- [Luke et al., Journal of Magnetic Resonance (2014)](http://dx.doi.org/10.1016/j.jmr.2014.01.011)
- [Spin Dynamics Wiki: `grad_pulse.m`](https://spindynamics.org/wiki/index.php?title=grad_pulse.m)
