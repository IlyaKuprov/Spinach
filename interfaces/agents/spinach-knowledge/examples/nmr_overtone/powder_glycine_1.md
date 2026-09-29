# examples/nmr_overtone/powder_glycine_1.m

- MATLAB implementation: [examples/nmr_overtone/powder_glycine_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/powder_glycine_1.m)

- Signature: `powder_glycine_1()`

## What the example models

This example calculates a powder-averaged 14N overtone NMR spectrum of glycine in the Fokker-Planck formalism, calling `powder(...,@overtone_pa,...,'qnmr')`. It sets `sys.magnet=14.1`, `sys.isotopes={'14N'}`, quadrupolar input `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`, and `inter.zeeman.scalar{1}=32.4`. Those are the source's literal parameter values; it gives no tensor units alongside them. The basis is full `sphten-liouv` (no approximation), with diagonal damping, zero equilibrium, and `damp_rate=500`.

The quadrupolar data are credited to O'Dell and Ratcliffe ([DOI 10.1016/j.cplett.2011.08.030](https://doi.org/10.1016/j.cplett.2011.08.030)). Unlike `mas_glycine_1.m`, this source does not define a rotor rate or MAS axis, rank limit, or phase multiplier: it uses the powder calculation route. It also explicitly cautions that the very short pulse uses unphysically large power. The source estimates calculation time as seconds. The citation is input-data provenance, not evidence that the plotted calculation reproduces an experiment.

## Powder average and pulse settings

The source uses grid `'rep_2ang_6400pts_sph'`, sweep `[0e3 15e3]`, 256 points and 256-point zero filling, with `axis_units='kHz'`. Preparation begins in 14N `Lz`; the coil and `Lx` operator are angle-weighted combinations of 14N `Lz` and `Lx`, with `theta=atan(sqrt(2))`. The pulse assignments are `rf_pwr=2*pi*11.3e6/sin(theta)`, `rf_dur=1e-6`, and `rf_frq=10e3`; the method is `'average'`.

## Output and limits

The script plots the real part of the spectrum with `plot_1d`. It does not save the spectrum, provide numerical peak positions or intensities, specify a contact time, or define a fitting workflow. Its pulse warning is part of the source and should be retained when interpreting this as an example configuration.
