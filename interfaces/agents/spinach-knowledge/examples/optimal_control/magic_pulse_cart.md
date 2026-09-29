# examples/optimal_control/magic_pulse_cart.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/magic_pulse_cart.m)

## Purpose

This Cartesian-control example frames a broadband 90° 13C excitation pulse as three simultaneous state transfers: Lz to Lx, Ly to Ly, and Lx to -Lz. Its model contains 100 non-interacting 13C spins distributed equally over chemical shifts from -100 to +100 ppm at 28.18 T. The accompanying reference is the magic-pulse paper at [doi:10.1016/j.jmr.2005.12.010](https://doi.org/10.1016/j.jmr.2005.12.010).

## Design objective and constraints

The example motivates broadband excitation that tolerates resonance offsets and RF power-calibration variation. The source gives a 50 microsecond duration ceiling from a worst-case 13C-1H coupling of about 200 Hz; the configured waveform has 40 one-microsecond intervals, or 40 microseconds total. It explores ten RF nutation levels between 50 and 70 kHz and optimises the Cartesian x/y controls with GRAPE, calling `fmaxnewton` with `@grape_xy` and `control.method='lbfgs'`. The initial control guess is constant in both quadratures. The configured iteration limit is 200; the source also names `NS` and `SNS` penalties with weights 0.01 and 10, respectively, without defining those labels in this example.

The three initial and target operators are normalised before optimisation. The basis configuration is `sphten-liouv` with `IK-2` at proximity level 1; the source comment describes this as retaining complete single-spin bases while omitting multi-spin orders for this case.

## Evaluation shown by the example

After optimisation, the waveform is simulated on an initial Lz state with the piecewise-constant exponential propagator. The resulting 13C signal is acquired with an L+ coil, a 70 kHz sweep, 2,048 points, and 16,384-point zero filling; a Gaussian apodisation parameter of 10 is applied before plotting the real spectrum on a ppm axis. A conventional hard-pulse spectrum is plotted for comparison; its configured power is `2*pi*60e3`, duration 4.2 microseconds, phase pi/2, and pulse rank 3.

These are design and evaluation settings in the example source, not a reported optimisation outcome or a kernel-test result. The source comments estimate a calculation time of minutes. Source contacts: ilya.kuprov@weizmann.ac.il and david.goodwin@inano.au.dk.
