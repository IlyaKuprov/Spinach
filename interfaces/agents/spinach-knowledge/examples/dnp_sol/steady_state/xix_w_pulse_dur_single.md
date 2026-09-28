# examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m

- Signature: `xix_w_pulse_dur_single()`

## Purpose

Calculates the steady-state proton signal for a XiX DNP experiment as a two-dimensional scan of electron-pulse duration and microwave offset. The source estimates a runtime of hours.

## Model and method

The W-band, 80 K model contains one trityl electron and one proton separated by 3.5 distance units. Proton relaxation is distance-dependent through `r1n_dnp`. The script forms the spin system in the untruncated Zeeman Liouville basis, detects proton (I_z), and evaluates the steady state with `xixdnp_steady` using an 800-point spherical powder grid. Each contact block is 360 ns and the shot spacing is 167 μs minus that contact time.

## Scan and output

The calculation covers 200 pulse durations from 2 to 21 ns and 101 microwave offsets from −230 to 205 MHz, at 20 MHz electron nutation frequency. The plot is updated after each duration and saved as `xix_w_pulse_dur_single.fig`.
