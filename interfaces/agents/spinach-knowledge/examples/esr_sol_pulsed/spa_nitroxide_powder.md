# examples/esr_sol_pulsed/spa_nitroxide_powder.m

- Function: `spa_nitroxide_powder()`.

## Model

This example calculates a soft-pulse powder EPR signal for a nitroxide radical. The source describes Fokker-Planck treatment of the soft pulse followed by time-domain acquisition and Fourier transform; its runtime estimate is seconds. Isotopes are `{'E','14N'}` and the electron Zeeman matrix is diagonal with values `[2.01045, 2.00641, 2.00211]`. The electron-`14N` coupling matrix is entered as `[1.2356 0 0.6322; 0 1.1266 0; 0.6322 0 8.2230]*1e7`. These interaction values and `sys.magnet=3.5` are reproduced as coded; the source does not annotate their units. The basis is `sphten-liouv` with no approximation, and trajectory-level SSR is disabled.

## Excitation, acquisition, and plotted result

The initial state and receiver operator are `Lz` and `L+` on electron `E`; `parameters.decouple={}`. The orientation grid is `rep_2ang_3200pts_sph`. Acquisition controls are offset `-2e8`, sweep `8e8`, 64 points, zero-fill 512, axis units MHz, derivative off, and axis inversion off. Raw offset and sweep units are not separately annotated.

The soft pulse is rank 2, phase `-pi/2`, frequency `-300e6`, duration `100e-9` s (100 ns), and power `2*pi*16.5e6`; the method is `expm`. The source gives no unit labels for the raw pulse-frequency and power values. It calls `powder(spin_system,@sp_acquire,parameters,'esr')`, applies `crisp` apodisation, computes `fftshift(fft(fid,parameters.zerofill))`, and plots the real part with `plot_1d`. No numerical spectrum is stated in the source.

## Source

[examples/esr_sol_pulsed/spa_nitroxide_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/spa_nitroxide_powder.m)
