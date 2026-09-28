# examples/dnp_sol/steady_state/tppm_q_rep_time_single.m

- Signature: `tppm_q_rep_time_single()`

## Purpose

Calculates the steady-state proton polarisation produced by a two-spin trityl–proton model under TPPM microwave irradiation, scanning the pulse repetition interval and plotting the detected proton `Lz` expectation value.

## Physical and numerical model

The model uses an electron and a proton at a 1.2142 T Q-band field, with trityl principal g values [2.00319, 2.00319, 2.00258], an 80 K spin temperature, and a 3.5 Å electron–nuclear separation. The relaxation model uses T1/T2 rates, including a distance- and orientation-dependent proton rate from `r1n_dnp`; the equilibrium is set to `dibari`, and only diagonal relaxation terms are retained. The spin system is represented in the full spherical-tensor Liouville formalism (`sphten-liouv`, no basis approximation).

## Experiment and output

The experiment detects proton `Lz`, uses an 800-point two-angle powder grid, and evaluates 30 logarithmically spaced repetition times from 10 μs to about 2 ms. Each point is computed by `powder(...,@xixdnp_steady,...,'esr')` with the source's TPPM pulse settings (including 120° second-pulse phase, −13 MHz added shift, and +2 MHz electron offset). The plotted real proton expectation value is saved as `tppm_q_rep_time_single.fig`.
