# examples/esr_sol_pulsed/eseem_nitroxide_powder.m

## Experiment and spin system

This example computes a powder-averaged two-pulse ESEEM signal for a `14N` nitroxide radical at 0.3249 T. The spin order is nitrogen then electron. The nitrogen self-coupling eigenvalues are `[-0.4, -1.6, +2.0] × 10^5 Hz`, with Euler angles `[0,0,0]`; the electron–nitrogen coupling is isotropic, `[2.0,2.0,2.0] × 10^6 Hz`, also with Euler angles `[0,0,0]`. The script does not explicitly assign an electron Zeeman tensor. The source comment names Figure 4a of [doi:10.1063/1.453532](https://doi.org/10.1063/1.453532) as the comparison target. Ideal hard pulses are assumed.

## Sequence and powder sampling

The calculation uses the `sphten-liouv` basis without approximation and disables trajectory-level SSR. `powder` calls the shared `eseem` sequence in the `esr` context using the finite `rep_2ang_400pts_sph` orientation grid. The initial state, receiver, screen, and pulse operator are electron `Lz`, `L+`, `L-`, and `Ly`, respectively; the pulse sequence is an ideal π/2 pulse, free evolution, ideal π pulse, and refocused evolution followed by receiver projection.

The sequence uses zero offset, 512 points, `timestep = 2e-7` s, and `zerofill = 2048`. Each of the two evolution intervals advances by `timestep/2`, hence the interpulse-delay increment is `1e-7` s and the full echo-time increment is `2e-7` s. These are fixed field, offset, and grid settings; the angular average is over the configured finite grid.

## Signal and displayed spectrum

Before the FFT, the script passes `mean(fid)-fid` to exponential apodisation with parameter 5. It zero-fills, Fourier-transforms, and applies `fftshift` to the result, then plots the real apodised time signal and the real spectrum. The spectral frequency axis uses `fft_freq_axis(npoints,timestep/2,zerofill-npoints)*1e-6` and is labelled in MHz. The script displays the figure; it contains no explicit data-file or figure-export call.

## Source links

- [Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/eseem_nitroxide_powder.m)
- [ESEEM sequence helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/eseem.m)
- [ESEEM helper reference](https://spindynamics.org/wiki/index.php?title=eseem.m)
- [Cited DOI](https://doi.org/10.1063/1.453532)
