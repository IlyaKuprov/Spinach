# examples/optimal_control/static_powder_control.m

Source: [examples/optimal_control/static_powder_control.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/static_powder_control.m)

- Signature: `static_powder_control()`

## Design objective and spin model

Design a deuterium pulse that prepares magnetisation in alanine's `-CD3` group for rephasing 100 microseconds after the pulse. The source describes a 600 MHz magnet, a single `2H` spin, and an alanine deuterium quadrupole-interaction tensor passed as `anas2mat(0,40e3,0,0,0,0)`. It uses the `sphten-liouv` basis without approximation. The normalised initial and target operators are deuterium `Lz` and `Lx`.

## Powder and robustness ensemble

The drift Liouvillians are built in the lab frame for the 100-orientation powder grid `rep_2ang_100pts_sph`. The settings specify no decoupled spins, zero central offset, and rotating frames `{{'2H',2}}`. RF power levels are `2*pi*[46,48,50,52,54]*1e3` rad/s (46–54 kHz when expressed as cycles per second), while five offset samples span −1 to +1 kHz. This tests the requested B1 and transmitter-offset spread across powder orientations.

## Pulse design

Goodwin's GRAPE Hessian method is selected, with 100 iterations, a 100 microsecond dead time, and `NS`/`SNS` penalties weighted 0.1 and 10. The pulse has 100 slices of 2 microseconds (200 microseconds total). A random 2-by-100 guess is optimised by `fmaxnewton` with `@grape_xy`; x/y components are then scaled by the mean power level.

## Observable and comparison

For every ensemble drift, the script applies the optimised pulse and then computes the target-observable evolution with 0.5 microsecond sampling and 499 intervals. It also computes an ideal free-induction reference initialised from the target state. The source plots the time-domain echo and compares Fourier transforms of the optimised half-echo (samples 201 onward) and the ideal FID (first 300 samples), with both spectra normalised by their own maxima. The header describes the design goal; no run output is included here to establish that rephasing was achieved. The source estimates minutes for calculation time.
