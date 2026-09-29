# examples/nmr_overtone/mas_valine_1.m

- MATLAB implementation: [examples/nmr_overtone/mas_valine_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/mas_valine_1.m)

- Signature: `mas_valine_1()`

## What the example models

This source models a 14N overtone MAS spectrum of N-acetylvaline with the Fokker-Planck formalism. It sets `sys.magnet=14.102` and `sys.isotopes={'14N'}`; the quadrupolar call is `eeqq2nqi(3.21e6,0.27,1,[0 0 0])`. It also supplies `inter.zeeman.eigs={[57.5 81.0 227.0]}` and `inter.zeeman.euler={[-90 -90 -17]*(pi/180)}`. These entries preserve the values and orientation expression supplied by the source; no tensor units are stated alongside them. The basis is `sphten-liouv` with no approximation. Relaxation uses diagonal damping, zero equilibrium, and `damp_rate=2000`.

The source attributes the valine quadrupolar tensor data to the authors' paper ([DOI 10.1039/C4CP03994G](https://doi.org/10.1039/C4CP03994G)) and estimates calculation time as hours. This attribution and estimate describe the example setup; they do not establish a fit or reproduction of an experimental spectrum.

## MAS settings and pulse-detected sequence

The source sets `theta=atan(sqrt(2))`, `max_rank=8`, axis `[sqrt(2/3) 0 sqrt(1/3)]`, rate `-19840`, and grid `'rep_2ang_6400pts_sph'`. It assigns sweep `[75e3 100e3]`, 256 points, 256-point zero filling, and `axis_units='kHz'`.

The initial state is 14N `Lz`; both the coil and the `Lx` operator are the angle-weighted `Lz/Lx` combinations using `cos(theta)` and `sin(theta)`. The RF settings are exactly `rf_pwr=2*pi*55e3/sin(theta)`, `rf_dur=70e-6`, and `rf_frq=86e3`; the method is `'average'`. The calculation calls `singlerot(...,@overtone_pa,...,'qnmr')`. The source gives no contact-time or fitting step.

## Output and distinction from the Z-detected variant

The result is multiplied by `exp(1i*1.75)` and `real(spectrum)` is plotted using `plot_1d`. The source does not save a spectrum file or report numerical output. This is the pulse-driven, angle-weighted preparation/detection path. In `mas_valine_2.m`, the stated initial state and coil are both 14N `Lz`, no RF pulse or phase multiplication is assigned, and the propagator is `@overtone_a`; those are distinct source-defined sequences, not interchangeable descriptions of one calculation.
