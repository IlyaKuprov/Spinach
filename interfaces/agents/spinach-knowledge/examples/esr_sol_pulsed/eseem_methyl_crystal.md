# examples/esr_sol_pulsed/eseem_methyl_crystal.m

## Experiment and spin system

This example computes a two-pulse X-band electron-spin-echo envelope modulation (ESEEM) response for a methyl radical at one fixed crystal orientation. The example imports vacuum-DFT magnetic parameters from `examples/standard_systems/methyl.log` with `g2spinach`; the Gaussian log contains four atoms (the carbon framework and three hydrogens), and the import maps the electron to `E` and the proton nuclei to `1H`. The imported model therefore carries the electron Zeeman and electron–proton hyperfine interactions from that calculation; their tensors are not literal assignments in this script. The field is 0.33 T. No powder average is performed. Ideal hard pulses are assumed.

## Sequence and sampling

The system is built in the `sphten-liouv` formalism with no basis approximation. The initial state is electron `Lz`, the receiver coil is electron `L+`, the screen is electron `L-`, and the pulse operator is electron `Ly`. `crystal` calls the shared `eseem` sequence in the `esr` context at orientation `[pi/5 pi/4 pi/3]` (radians). The helper applies an ideal π/2 pulse, evolves for an interpulse interval, applies an ideal π pulse, evolves for the refocused interval, and projects onto the receiver.

The example requests 512 points, `timestep = 1e-8` s, and zero-fills to 4096 points. In the ESEEM helper each of the two evolution periods advances by `timestep/2` per sample, so the interpulse-delay increment is 5 ns and the full echo-time increment is 10 ns. The sequence is a single-orientation calculation at a fixed field and orientation, not a parameter sweep.

## Signal and displayed spectrum

The upper panel plots the real FID against sample index times `timestep`, labelled in microseconds. For the spectral panel the script subtracts the FID mean, applies Kaiser apodisation with parameter 6, computes a 4096-point FFT, applies `fftshift`, and plots its magnitude. Its frequency axis is constructed directly as `linspace(-1/timestep,1/timestep,zerofill)*1e-6` and labelled in MHz; this is the axis expression used by the example. The script plots the panels but contains no explicit data-file or figure-export call.

## Source links

- [Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/eseem_methyl_crystal.m)
- [Imported methyl DFT log](https://github.com/IlyaKuprov/Spinach/blob/main/examples/standard_systems/methyl.log)
- [ESEEM sequence helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/eseem.m)
- [ESEEM helper reference](https://spindynamics.org/wiki/index.php?title=eseem.m)
