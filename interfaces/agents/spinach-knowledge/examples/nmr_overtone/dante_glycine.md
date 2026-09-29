# examples/nmr_overtone/dante_glycine.m

- MATLAB implementation: [examples/nmr_overtone/dante_glycine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/dante_glycine.m)

## Purpose and literature provenance

This example calculates a 14N overtone DANTE spectrum of glycine with the Fokker-Planck formalism. The source attributes the glycine quadrupolar tensor data to O'Dell and Ratcliffe, https://doi.org/10.1016/j.cplett.2011.08.030, and estimates minutes of calculation time. The paper is cited as provenance for the tensor input; the spectrum and its phase adjustment below are simulation output, not a claim of reproducing a measured spectrum.

## Spin model and MAS

The source defines a single 14N spin with sys.magnet=14.1. Its coupling input is eeqq2nqi(1.18e6,0.53,1,[0 0 0]), and its scalar Zeeman input is 32.4. These are transcribed without assigning units beyond those stated by the source. Relaxation is damp, retaining diagonal terms with zero equilibrium and damp_rate=300. The basis is sphten-liouv with no approximation; krylov and trajlevel are disabled.

The magic-angle parameter is atan(sqrt(2)); the spectrum axis input is [sqrt(2/3) 0 sqrt(1/3)]. The simulation sets max_rank=7, rotor-rate input -19840, and grid rep_2ang_1600pts_sph. The sweep is [-60e3 80e3] with 2048 points and 2048-point zero-fill, and axis_units=kHz. The initial state is the 14N Lz state; receiver and RF operator use magic-angle-weighted Lz and Lx terms.

## DANTE settings and output

The sequence uses average treatment, four periods, and two pulses. The source sets pulse amplitude to 2*pi*55e3/sin(theta), pulse duration to 10e-6, and RF-frequency input to 48e3; it selects 14N as the irradiated spin. The simulated spectrum is calculated with singlerot and overtone_dante, then multiplied by exp(-1i*2.12) before its real part is plotted with plot_1d.

This is a single-spin nitrogen overtone DANTE calculation, not a proton-to-nitrogen cross-polarisation scan. The phase multiplication is part of the source's plotting workflow; it is not evidence of a fitted experimental phase.
