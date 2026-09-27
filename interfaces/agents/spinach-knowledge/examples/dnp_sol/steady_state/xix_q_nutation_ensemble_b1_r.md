# examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m

- Signature: `xix_q_nutation_ensemble_b1_r()`

## Purpose

Compares steady-state XiX DNP field profiles at six electron nutation frequencies while averaging over electron–proton distance and microwave B1 distributions. The source estimates a calculation time of minutes.

## Physical and numerical setup

The six nutation frequencies are 6.8, 9.6, 13.5, 17.5, 25, and 36 MHz, each paired with its listed shot repetition time (0.051, 0.051, 0.102, 0.153, 0.153, and 0.306 ms). For each frequency, the script samples distance from 3.5–20 Å with three Gauss–Legendre nodes and B1 from 0.2 to 1.2 times that frequency with five nodes. The spin system is an electron–proton pair at 80 K and 1.2142 T; the XiX pulse train has 36 blocks and 48 ns pulses.

## Calculation and output

For each distance/B1 pair, `powder(...,@xixdnp_steady,...,'esr')` evaluates the steady state at 13 offsets from −64 to −52 MHz. The results are averaged over B1 weights and over the radial distance distribution, including the (r^2) radial Jacobian. The script plots the negative real proton (I_z) expectation versus offset and nutation frequency, then saves `xix_q_nutation_ensemble_b1_r.fig`.
