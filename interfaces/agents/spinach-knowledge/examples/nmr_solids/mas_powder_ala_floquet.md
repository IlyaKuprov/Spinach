# examples/nmr_solids/mas_powder_ala_floquet.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_ala_floquet.m

## What it models

The example computes a 13C MAS powder spectrum for alanine. Its header says to assume 1H decoupling and identifies the Floquet MAS formalism. The spin-system input is parsed from ../standard_systems/alanine.log, described in the source as a PCM DFT calculation, and passed through g2spinach with the C-to-13C and N-to-15N mapping and the numeric pair [182.1 264.5]. The field assignment is sys.magnet=14.1. The source does not state units alongside the two mapping values.

The basis is sphten-liouv with approximation none, longitudinal 15N, and projection +1. Inter- and proximity-cutoff values are set to 5.0 and 4.0, respectively. The acquisition configuration sets rate=2000, axis=[1 1 1], max_rank=17, sweep=5e4, npoints=256, zerofill=1024, offset=15000, and the 13C spin. It also sets decouple={} and a rep_2ang_100pts_sph grid. It sets axis_units=ppm and invert_axis=1. The source does not attach units to rate, sweep, offset, or cutoff values. In particular, its header's 1H-decoupling assumption is not accompanied here by a proton decoupling pulse sequence; the empty decouple assignment is reported as written, not interpreted as evidence for how decoupling is implemented.

## Calculation and displayed signal

The initial state and receiver are both 13C L+. The example calls floquet with @acquire in nmr mode, applies exponential apodisation with parameter 6, Fourier transforms the FID with the 1024-point zero-fill, shifts the spectrum, and plots its real part with plot_1d. The source comment estimates a calculation time of minutes; it does not report a measured spectrum or a separate validation result. No DOI is supplied in the source or the inspected page baseline.
