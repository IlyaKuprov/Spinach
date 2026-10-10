# examples/optimal_control/features_bss.m

[Stable source link](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_bss.m) · [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=features_bss.m)

## Calculation

This is a simulated 90° proton-pulse design in the regime where Bloch–Siegert (counter-rotating-field) effects matter. The one-spin model is on resonance at 1 MHz: `sys.magnet=2π×10⁶/spin('1H')` sets that Larmor frequency and the chemical shift is 0 ppm. Spinach builds the full single-spin `sphten-liouv` basis with no approximation. The normalised initial and target states are longitudinal magnetisation (`Lz`) and transverse magnetisation (`Lx`). The drift is set to zero; the two quadrature controls are `Lx` and `Ly` for 1H.

The ensemble varies an offset through the transverse `Lx` operator, using the 11 values in `linspace(-1e5,1e5,11)`. The source does not attach units to these values, so they are reported as entered rather than relabelled as Hz. Its RF scale is 0.2 times the absolute proton base frequency, and it uses 50 equal slices with `pulse_dt=(8π/pwr_levels/50)`; using the source's 1 MHz base frequency gives a total duration of 20 μs. The deterministic initial waveform has 50 amplitude samples ramping from 0.1 to 0.5 and 50 second-quadrature samples at 0.05.

## Optimisation and reported observable

`optimcon` prepares the control problem for LBFGS GRAPE (`fmaxnewton` with `@grape_xy`), real-valued fidelity and a 100-iteration limit. The script optimises once with `control.bsiegert=true` and again with it false, then evaluates each resulting pulse using `ensemble` in the corrected model. It reports two corrected-model ensemble fidelities without hard-coded values. The source header describes the corrected-versus-uncorrected gap as about six times the Lz-offset ensemble gap; the script does not compute that separate offset comparison. The calculation is a Spinach simulation, not a hardware measurement.

The zero drift uses the compiled state-space dimension `bas.offsets(end)`, matching the state and control-operator dimensions.
