# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1_r.m

- Signature: `xix_q_con_time_ensemble_b1_r()`

## Purpose

Calculates steady-state proton polarisation versus XiX contact time, averaging over both electron–proton distance and electron microwave-field ensembles.

## Model and quadrature

The system is a trityl electron and proton at 1.2142 T and 80 K. Three Gauss–Legendre distance nodes span 3.5–20 Å; five microwave-field nodes span 10–20 MHz. At each distance, the source sets the pair coordinates and updates the orientation-dependent proton T1 rate through `r1n_dnp`; T2 rates, diagonal relaxation retention and `dibari` equilibrium are specified. The full spherical-tensor Liouville basis is used without basis approximation.

## XiX scan and output

For each distance, field node and loop count (1–64), the calculation uses a 48 ns pulse, inverted second-pulse phase, an 800-point two-angle spherical powder grid, and `powder(...,@xixdnp_steady,...,'esr')`. Shot spacing is 153 μs minus the total pulse duration; the source sets a −13 MHz added shift and +61 MHz electron offset. The signal is weighted over microwave-field nodes and over distance nodes with the radial `r^2` Jacobian. The resulting real proton `Lz` expectation value is plotted against total contact time and saved as `xix_q_con_time_ensemble_b1_r.fig`.
