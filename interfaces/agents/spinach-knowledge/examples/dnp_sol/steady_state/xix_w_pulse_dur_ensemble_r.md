# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_r.m

- Signature: `xix_w_pulse_dur_ensemble_r()`

## Purpose

Computes the steady-state proton signal for a XiX DNP experiment over electron-pulse duration and microwave offset, averaged over an electron–proton distance ensemble. The source estimates a runtime of hours.

## Model and method

The W-band model uses one electron and one proton at 80 K, with the trityl electron g tensor and a proton shift. Three Gauss–Legendre distance nodes are generated over 3.5–20 (in the source's distance units); at each node the proton relaxation rate is evaluated using `r1n_dnp`. A spherical powder average calls `xixdnp_steady` with a 20 MHz electron nutation frequency, 800-point grid, and a 360 ns contact block.

## Scan and output

For each distance, the script calculates 200 pulse durations from 2 to 21 ns and 101 microwave offsets from −230 to 205 MHz. It combines the resulting proton (I_z) signals with the quadrature weights and the (r^2) radial Jacobian, then plots signal versus pulse duration and offset. The figure is saved as `xix_w_pulse_dur_ensemble_r.fig`.
