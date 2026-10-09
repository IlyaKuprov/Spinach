# examples/nmr_overtone/mas_valine_2.m

- MATLAB implementation: [examples/nmr_overtone/mas_valine_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/mas_valine_2.m)

- Signature: `mas_valine_2()`

## What the example models

This is the Z-detected 14N overtone MAS example for N-acetylvaline, evaluated with the Fokker-Planck formalism and `singlerot(...,@overtone_a,...,'qnmr')`. Its spin and interaction model matches the values explicitly assigned in `mas_valine_1.m`: `sys.magnet=14.102`, `sys.isotopes={'14N'}`, quadrupolar input `eeqq2nqi(3.21e6,0.27,1,[0 0 0])`, Zeeman eigenvalues `[57.5 81.0 227.0]`, and Euler expression `[-90 -90 -17]*(pi/180)`. The source does not state tensor units next to these values. The full `sphten-liouv` basis has no approximation; relaxation is diagonal damping with zero equilibrium and `damp_rate=2000`.

The source attributes the valine quadrupolar tensor data to the authors' paper ([DOI 10.1039/C4CP03994G](https://doi.org/10.1039/C4CP03994G)) and estimates calculation time as hours. These comments do not report a fit, experimental reconstruction, or numerical validation.

## MAS and Z-detection

The source assigns `max_rank=8`, axis `[sqrt(2/3) 0 sqrt(1/3)]`, rotor-rate parameter `-19840`, and grid `'rep_2ang_6400pts_sph'`. The sweep is `[75e3 100e3]`, with 256 points, 256-point zero filling, and `axis_units='kHz'`.

Both the initial state and the coil are the 14N `Lz` state. The source sets `spins={'14N'}` and does not assign RF-pulse fields, an angle-weighted `Lx` detector, or a phase factor. It invokes `singlerot` with `@overtone_a` and `'qnmr'`; the method is therefore distinct from the explicitly pulse-driven `mas_valine_1.m` path. No contact-time or fitting procedure is specified.

## Output and limits

The script plots `real(spectrum)` with `plot_1d`; it does not save a spectrum file or state its numerical values. The absence of pulse and phase assignments here is a source-level distinction, not a claim about how a measured experiment was acquired.
