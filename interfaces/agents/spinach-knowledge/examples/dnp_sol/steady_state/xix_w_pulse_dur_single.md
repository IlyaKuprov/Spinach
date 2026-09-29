# examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m)

## Purpose

Runs a steady-state XiX dynamic nuclear polarisation (DNP) scan over electron-pulse duration and microwave resonance offset for an electron–proton pair. The source estimates a calculation time of hours.

## Run and model

Call `xix_w_pulse_dur_single()` from MATLAB. The calculation builds a W-band model with one electron and one proton, then calls `powder(spin_system,@xixdnp_steady,parameters,'esr')` at each pulse duration. The reported model uses the trityl electron Zeeman values `[2.00319 2.00319 2.00258]`, proton values `[0 0 5]` (the source describes the trityl values as a g-tensor and the 1H values as ppm guesses), Euler angles `[0 10 0]` and `[0 0 10]` degrees, spin temperature 80, and coordinates `[0 0 0]` and `[0 0 3.5]`. The source does not state coordinate units. The field setting is `sys.magnet=3.4` (described as W-band).

The basis is `sphten-liouv` with no approximation; propagator chopping tolerance is `1e-12`. Relaxation uses `r1n_dnp`, an electron rate of `1e3`, `r2_rates={200e3 50e3}`, diagonal relaxation retention, and `dibari` equilibrium. The electron–proton coordinate separation is passed to the R1 model as `r_en`.

## Dependencies

Uses Spinach `create`, `basis`, `state`, and `powder`, the `xixdnp_steady` sequence callback, and `r1n_dnp` for the distance-dependent R1 rate. Plotting and figure output use `kfigure`, `imagesc`, the Spinach axis/colorbar helpers, and MATLAB `savefig`.

## Scan and sequence

Electron pulse duration is sampled at 200 points from 2 to 21 ns. The 101-point offset grid runs from `-230e6` to `205e6` (the plotted axis is MHz); `addshift=-33e6`. Electron nutation frequency is `20e6` Hz, and the spherical orientation grid is `rep_2ang_800pts_sph`. For each duration `τ`, the sequence uses `round(360e-9/(2τ))` loops and a shot spacing of `167e-6 - 360e-9`. Its distinctive XiX setting is `phase=pi`, annotated in the source as the second pulse having inverted phase.

## Output and scope

The powder calculation fills a complex-valued array over offset and pulse duration; the image displays its real part as the proton `I_z` expectation value. The axes are microwave resonance offset (MHz) and pulse duration (ns). The figure is saved as `xix_w_pulse_dur_single.fig` in the current working directory. This is the source's finite two-parameter grid and steady-state XiX callback, not a result for unscanned settings.
