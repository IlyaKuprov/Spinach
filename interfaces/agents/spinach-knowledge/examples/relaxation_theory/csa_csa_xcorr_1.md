# examples/relaxation_theory/csa_csa_xcorr_1.m

- Signature: `csa_csa_xcorr_1()`
- Source: [examples/relaxation_theory/csa_csa_xcorr_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/csa_csa_xcorr_1.m)

This example constructs and displays a Bloch–Redfield–Wangsness (BRW) relaxation superoperator for a two-spin `1H`–`13C` system. Its two anisotropic shielding interactions generate CSA relaxation, including CSA–CSA cross-correlation; the source notes that Spinach's relaxation module includes all cross-correlations for this model.

The field is `14.1 T`. The shielding principal values are `[7 15 -22]` and `[11 18 -29]` ppm, with Euler orientations `[pi/5 pi/3 pi/11]` and `[pi/6 pi/7 pi/15]` radians. The model selects Redfield relaxation, a single correlation time `tau_c=2e-9 s`, zero equilibrium, and `labframe` retention. The Liouville-space basis is `sphten-liouv` with `approximation='none'`.

A single correlation time specifies a one-timescale tumbling model, not an anisotropic diffusion tensor. In the usual isotropic rotational-diffusion interpretation, rank-2 interaction correlations decay exponentially and give a Lorentzian spectral density proportional to `tau_c/[1 + (omega*tau_c)^2]`; this is the frequency weighting used by BRW relaxation. The example does not select individual CSA cross terms manually: it asks `relaxation(spin_system)` to assemble them from the specified interactions.

The function prints the full relaxation matrix with `disp(full(...))`. It does not propagate an initial state, define detection, generate a spectrum, or compare with experimental data; the output is a superoperator, not a measured or simulated line shape.
