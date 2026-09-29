# examples/esr_sol_pulsed/eseem_phenyl_crystal.m

## Experiment and spin system

This example calculates the two-pulse X-band ESEEM response of a phenyl radical at a single crystal orientation. It imports the vacuum-DFT magnetic parameters from `examples/standard_systems/phenyl.log` through `g2spinach`, mapping the electron to `E` and proton nuclei to `1H`. The Gaussian log contains eleven atoms, corresponding to the phenyl framework and its five hydrogens; the imported interaction tensors are not entered as literal values in this example. The script sets the field to 0.33 T and assumes ideal hard pulses.

## Sequence and sampling

The spin system uses the `sphten-liouv` basis without approximation. Electron `Lz` is the initial state, electron `L+` and `L-` are the receiver and screen, and electron `Ly` is the pulse operator. `crystal` calls the shared `eseem` sequence in the `esr` context with fixed orientation `[pi/5 pi/4 pi/3]` (radians); this is not a powder average. The sequence applies an ideal π/2 pulse, evolves for an interpulse interval, applies an ideal π pulse, and evolves through the refocused interval before receiver projection.

The run requests 512 points with `timestep = 1e-8` s and zero-fills to 4096. The helper advances each of the two evolution periods by `timestep/2`, so the interpulse-delay increment is 5 ns and the full echo-time increment is 10 ns. The source comments give a calculation time of minutes.

## Signal and displayed spectrum

The upper panel plots the real FID against sample index times `timestep` in microseconds. The script then removes the FID mean, applies Kaiser apodisation with parameter 6, computes a 4096-point FFT, applies `fftshift`, and plots its magnitude. The frequency axis is `fft_freq_axis(npoints,timestep/2,zerofill-npoints)*1e-6`, labelled in MHz. The script displays the two panels and contains no explicit data-file or figure-export call.

## Source links

- [Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/eseem_phenyl_crystal.m)
- [Imported phenyl DFT log](https://github.com/IlyaKuprov/Spinach/blob/main/examples/standard_systems/phenyl.log)
- [ESEEM sequence helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/eseem.m)
- [ESEEM helper reference](https://spindynamics.org/wiki/index.php?title=eseem.m)
