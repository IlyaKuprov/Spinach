# examples/esr_sol_pulsed/eseem_nitroxide_crystal.m

## Experiment and spin system

This is a two-pulse X-band ESEEM calculation for a nitroxide radical at one specified crystal orientation. The script explicitly constructs two spins, `E` and `14N`, at 0.33 T. Its electron Zeeman tensor is diagonal with principal values `gxx = 2.01045`, `gyy = 2.00641`, and `gzz = 2.00211`. The electron–nitrogen coupling matrix is entered in the source as `1e7` times `[[1.2356, 0, 0.6322], [0, 1.1266, 0], [0.6322, 0, 8.2230]]` (Hz in Spinach frequency units). No nitrogen self-coupling tensor is assigned in this example. The source describes the magnetic parameters as DFT-derived, but gives these matrix values directly rather than importing them from a log. Ideal pulses are assumed.

## Sequence and sampling

The basis is `sphten-liouv` without approximation. The initial state is electron `Lz`; electron `L+` and `L-` provide the coil and screen, and electron `Ly` is the pulse operator. `crystal` calls the shared `eseem` sequence in the `esr` context at fixed orientation `[pi/5 pi/4 pi/3]` (radians); there is no powder average.

The run uses 1024 points, `timestep = 1.25e-8` s, and `zerofill = 4096`. The helper evolves each of the two echo periods by `timestep/2`, giving a 6.25 ns interpulse-delay increment and a 12.5 ns full echo-time increment. It applies an ideal π/2 pulse, the first free-evolution interval, an ideal π pulse, and the refocused interval before detection.

## Signal and displayed spectrum

The upper panel plots the real time-domain signal against sample index times `timestep` in microseconds. Before the spectral transform the script subtracts the signal mean and applies Kaiser apodisation with parameter 6. It then computes a 4096-point FFT, applies `fftshift`, and plots its magnitude. The frequency axis uses `fft_freq_axis(npoints,timestep/2,zerofill-npoints)*1e-6`, labelled in MHz, so the axis is based on the interpulse-delay increment. The script displays the figure and has no explicit data-file or figure-export call.

## Source links

- [Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/eseem_nitroxide_crystal.m)
- [ESEEM sequence helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/eseem.m)
- [ESEEM helper reference](https://spindynamics.org/wiki/index.php?title=eseem.m)
