# examples/relaxation_theory/csa_csa_xcorr_2.m

- Signature: `csa_csa_xcorr_2()`
- Source: [examples/relaxation_theory/csa_csa_xcorr_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/csa_csa_xcorr_2.m)

This liquid-state NMR example shows how CSA–CSA cross-correlation in a pair of `103Rh` shielding tensors affects the widths of the three lines in a proton triplet. The source describes the output in those terms; it is a Spinach simulation, not evidence of agreement with a particular experiment.

The spin system is `1H, 103Rh, 103Rh` at `11.75 T`. The proton shielding matrix is isotropic, with diagonal values `[6.9 6.9 6.9]` ppm; both rhodium matrices are `diag([7250 8000 7250])` ppm. Scalar couplings are `4 Hz` from the proton to each rhodium and `100 Hz` between the rhodium nuclei. Redfield relaxation uses one correlation time, `tau_c=10e-9 s`, zero equilibrium, and `rlx_keep='secular'`, retaining the secular relaxation terms. The basis is `sphten-liouv` without approximation. A single `tau_c` represents a one-timescale rotational-correlation model; under isotropic tumbling the rank-2 correlation spectrum is Lorentzian, with frequency dependence proportional to `tau_c/[1 + (omega*tau_c)^2]`.

For acquisition, the source sets both the initial state and receiver to the proton `L+` operator, with no decoupling. It specifies offset `6.9*500 Hz`, sweep width `50 Hz`, `2048` points, zero-fill to `16384`, ppm axis units, and an inverted axis. `liquid(...,@acquire,...,'nmr')` generates the FID; exponential apodisation with coefficient `20` precedes the Fourier transform and plot. These are simulation settings, not experimental acquisition metadata.
