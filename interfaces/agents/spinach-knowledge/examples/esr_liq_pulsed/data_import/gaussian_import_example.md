# examples/esr_liq_pulsed/data_import/gaussian_import_example.m

- MATLAB implementation: [examples/esr_liq_pulsed/data_import/gaussian_import_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/data_import/gaussian_import_example.m)

## Purpose

`gaussian_import_example()` imports Gaussian output for a methyl-radical liquid ESR simulation. Its source comment associates the unusual signal-intensity pattern with cross-correlation between the g tensor and hyperfine couplings.

## Input and imported system

The function expects `gaussian_methyl_radical.out` to be accessible to `gparse`. It passes the parsed output to `g2spinach` with isotope mapping `{{'E','E'},{'H','1H'}}`, the `[0 0]` argument, and `options.no_xyz=1`. It then sets `sys.magnet=0.339`; the source does not state a unit for this setting.

## ESR model and acquisition

The basis is `sphten-liouv` with `approximation='none'`. Relaxation is configured as Redfield, equilibrium as `'zero'`, retained relaxation terms as `'secular'`, and `inter.tau_c={5e-10}`. The source provides no unit next to that correlation-time setting. Initial and detected states are electron `L+`; the sequence parameters select `{'E'}` and leave `decouple` empty.

The exact simulation call is `liquid(spin_system,@acquire,parameters,'esr')`. Acquisition settings are `sweep=5e8`, `npoints=512`, `zerofill=1024`, `axis_units='mT'`, `derivative=1`, and `invert_axis=0`. The script does not label the numeric sweep value with a unit, so it is best retained as the configured value rather than restated as a bandwidth in a different unit.

## Processing, dependencies, and scope

The FID is passed through `apodisation` with `{{'none'}}`, transformed as `fftshift(fft(fid,parameters.zerofill))`, and plotted as the real spectrum. With Spinach on the MATLAB path and the named Gaussian output available, call `gaussian_import_example()`. The import path uses Spinach's `gparse`/`g2spinach` interfaces; simulation and display use `create`, `basis`, `state`, `liquid`, `acquire`, `apodisation`, and the Spinach plotting helpers. It opens a plot and declares no return value; the FID and spectrum are local variables. The example demonstrates this specific imported parameter set and acquisition configuration, not a general Gaussian-output compatibility guarantee.
