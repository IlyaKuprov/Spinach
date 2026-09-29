# examples/nmr_overtone/mas_glycine_1.m

- MATLAB implementation: [examples/nmr_overtone/mas_glycine_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/mas_glycine_1.m)

- Signature: `mas_glycine_1()`

## What the example models

This is a single-spin 14N overtone MAS calculation using the Fokker-Planck formalism and `singlerot(...,@overtone_pa,...,'qnmr')`. The source sets `sys.magnet=14.1`, `sys.isotopes={'14N'}`, and constructs the quadrupolar interaction with `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`; it sets `inter.zeeman.scalar{1}=32.4`. These are the source's literal values and constructor arguments; the source does not attach additional units to the tensor or scalar values. The basis is the full `sphten-liouv` basis (`approximation='none'`). Relaxation is `{'damp'}`, with diagonal terms retained, zero equilibrium, and `damp_rate=300`.

The source cites O'Dell and Ratcliffe for the glycine quadrupolar tensor data ([DOI 10.1016/j.cplett.2011.08.030](https://doi.org/10.1016/j.cplett.2011.08.030)). Separately, it says its simulation parameters are set to reproduce Figure 3b of the authors' paper ([DOI 10.1039/C4CP03994G](https://doi.org/10.1039/C4CP03994G)). That is a statement of the example's intended comparison, not evidence here of an experimental reproduction or an independently checked match.

## Rotor, preparation and detection

The source defines the magic angle as `atan(sqrt(2))`. The single-rotor settings are rank 6, axis `[sqrt(2/3) 0 sqrt(1/3)]`, rate `-19840`, and grid `'rep_2ang_6400pts_sph'`. The sweep is `[44e3 52e3]`, with `axis_units='kHz'`, 256 points and 256-point zero filling; these are the script's assigned values, not converted values.

Preparation starts from `state(...,'Lz','14N')`. Both the coil and the `Lx` operator are the same angle-weighted combination of the 14N `Lz` and `Lx` states/operators, with coefficients `cos(theta)` and `sin(theta)`. The pulse fields are `rf_pwr=2*pi*55e3/sin(theta)`, `rf_dur=260e-6`, and `rf_frq=48e3`; the method is `'average'`. No contact-time or fitting workflow is defined in this source.

## Output and limits

The script phases the returned spectrum by `exp(1i*1.35)` and plots `real(spectrum)` with `plot_1d`. The source labels the calculation time as minutes. It defines a plotted simulation, not a saved data product; it does not report a fit, validation, or the numerical spectrum values.
