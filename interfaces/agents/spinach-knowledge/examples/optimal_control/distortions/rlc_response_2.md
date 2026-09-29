# examples/optimal_control/distortions/rlc_response_2.m

Source: [examples/optimal_control/distortions/rlc_response_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/rlc_response_2.m)

- Signature: `rlc_response_2()`
- Calculation time: the source header estimates minutes.

## Purpose and model

This example asks how a probe-circuit response model changes a simulated deuterium pre-phasing pulse for the CD3 group of alanine. The source describes the goal as placing the deuterium magnetisation for rephasing 100 microseconds after the pulse. Its spin model contains 2H, with the alanine-CD3 NQI input `anas2mat(0,40e3,0,0,0,0)`; the script labels the magnet as 600 MHz and sets `sys.magnet=14.1` T. It uses the full sphten-liouv basis (`bas.approximation='none'`) and a 100-orientation powder grid named `rep_2ang_100pts_sph`. The powder drift ensemble is built with `drifts(...,@powder,...,'labframe')`.

## Pulse design and response

The normalised initial and target states are the 2H `Lz` and `Lx` states. The two controls are the 2H `Lx`/`Ly` operators on one channel. The amplitude-robust design spans five RF levels, 40, 45, 50, 55, and 60 kHz per channel (stored as `2*pi*[40 45 50 55 60]*1e3`). It uses 75 slices of 2 microseconds each, a 100-microsecond optimiser dead time, a maximum of 50 iterations, an NS penalty of weight 1, and freezes the first and last four slices. The source sets the final optimiser method to `goodwin`, with the rectangle integrator, and calls `grape_xy` from `randn(2,75)/10` after setting the first and last four samples to zero. Optimiser plots are requested for XY controls, robustness, and the spectrogram.

The resulting Cartesian control arrays are scaled by the mean power level and passed to `restrans` with the 2H Larmor frequency at 14.1 T, Q = 200, and the `pwc` response option (the call also supplies a final argument of 100). This produces the circuit-response figure; the page describes a model calculation, not a measured probe trace. The script contains no reported response values or comparison metric, so it does not establish hardware performance or a numerical loss of pulse accuracy.
