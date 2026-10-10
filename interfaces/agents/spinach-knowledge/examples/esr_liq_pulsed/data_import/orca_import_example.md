# examples/esr_liq_pulsed/data_import/orca_import_example.m

- MATLAB implementation: [examples/esr_liq_pulsed/data_import/orca_import_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/data_import/orca_import_example.m)

## Purpose

`orca_import_example()` imports ORCA output for a methyl-radical liquid ESR simulation. Its source comment associates the unusual signal-intensity pattern with cross-correlation between the g tensor and hyperfine couplings.

## Input and imported system

The function expects `orca_methyl_radical.out` to be accessible to `oparse`. It reads the ORCA properties, defines one electron and three protons (`{'E','1H','1H','1H'}`), and assigns the electron g tensor from `props.g_tensor.matrix`. The three electron–proton coupling entries are populated from `props.hfc.full.matrix{2}`, `{3}`, and `{4}` using the source expression `1e6*gauss2mhz(...)`. This is the literal conversion used by the example; the page does not reinterpret its resulting unit. The magnet setting is `sys.magnet=0.339`, with no unit stated in the source.

## ESR model and acquisition

The basis is `sphten-liouv` with `approximation='none'`. Relaxation is configured as Redfield, equilibrium as `'zero'`, retained relaxation terms as `'secular'`, and `inter.tau_c={5e-10}`; the script gives no unit alongside that value. Initial and detected states are electron `L+`; the sequence parameters select `{'E'}` and leave `decouple` empty.

The simulation call is `liquid(spin_system,@acquire,parameters,'esr')`. Acquisition settings are `sweep=5e8`, `npoints=512`, `zerofill=1024`, `axis_units='mT'`, `derivative=1`, and `invert_axis=0`. The numeric sweep value has no unit comment in this source.

## Processing, dependencies, and scope

The FID uses `apodisation` with `{{'none'}}`, then `fftshift(fft(fid,parameters.zerofill))`; the plotted spectrum is its real part. With Spinach on the MATLAB path and the named ORCA output available, call `orca_import_example()`. The import uses `oparse` and `gauss2mhz`; simulation and display use Spinach system/basis/state construction, `liquid`, `acquire`, `apodisation`, and plotting helpers. The function opens a plot and declares no return value; FID and spectrum are local. This page describes the example's explicit tensor-transfer route and configured input, not universal coverage of ORCA output variants.
