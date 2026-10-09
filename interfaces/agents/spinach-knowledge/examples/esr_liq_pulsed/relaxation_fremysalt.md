# examples/esr_liq_pulsed/relaxation_fremysalt.m

## Purpose and interface

A pulse-acquire FFT ESR version of the EasySpin Fremy-salt test example, with acknowledgement to Stefan Stoll. It uses explicit Liouville-space time propagation and Redfield relaxation and is set to reproduce Figure 3a of the cited paper. Call `relaxation_fremysalt()` with no arguments; the function declares no return values and displays a plot rather than saving data.

## Spin system and relaxation

The spin system is one electron and one `14N` nucleus. The electron g principal values are `[2.00785 2.00590 2.00265]` with zero Euler angles. The electron-nitrogen coupling principal values are `[15.4137 14.0125 80.4316]*1e6`, also with zero Euler angles; the source does not label coupling units. The magnetic-field setting is `sys.magnet=0.33` (no unit is written beside the assignment). The basis is full `sphten-liouv`. Relaxation is Redfield with zero equilibrium, secular terms retained, and `inter.tau_c={8e-10}` (the source gives no unit for this value).

## Pulse-acquire and processing

The initial state and receiver are the electron `L+` operator; decoupling is empty and offset is `-2e7`. Acquisition uses sweep `2e8`, 512 points, 1024-point zero filling, GHz-labframe axis units, derivative display, and inverted axis. It calls `liquid(...,@acquire,...,'esr')`, applies no apodisation, Fourier transforms the FID, and plots the real spectrum. No numerical spectrum is embedded in the source; the figure is the runtime output.

Requires the Spinach system/basis, ESR acquisition, processing, and plotting functions (`create`, `basis`, `state`, `liquid`, `acquire`, `apodisation`, `kfigure`, `plot_1d`). Reference: [Figure 3a, DOI: 10.1209/epl/i2004-10459-y](https://doi.org/10.1209/epl/i2004-10459-y).

[Source: `relaxation_fremysalt.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/relaxation_fremysalt.m).
