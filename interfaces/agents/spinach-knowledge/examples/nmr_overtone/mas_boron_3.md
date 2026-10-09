# examples/nmr_overtone/mas_boron_3.m

- MATLAB implementation: [examples/nmr_overtone/mas_boron_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/mas_boron_3.m)

This example sets up a panoramic 10B MAS overtone spectrum. The source comments describe JEOL-direction spinning, attribute parameter choices to Nghia Duong and Yusuke Nishiyama, and say an unphysically strong pulse is used to obtain a panoramic spectrum; the comments estimate hours of calculation. No DOI, experimental result, or fit comparison is given in the source.

The model uses isotope 10B, magnet setting 16.4, and coupling input `eeqq2nqi(0.7e6,0.0,3,[0 0 0])`. The source does not annotate units for the magnet or coupling arguments. Relaxation is `damp` with diagonal retention, zero equilibrium, and `damp_rate=1000`. The basis is `sphten-liouv` with approximation `none`.

The sequence parameters are rank 12, axis `[sqrt(2/3) 0 sqrt(1/3)]`, rate 70000, and grid `rep_2ang_200pts_oct`. The spectrum uses sweep `[-200e3 200e3]`, 4096 points, 4096-point zero-fill, and `axis_units='kHz'`. Both the initial state and receiver are 10B `Lz`. The code calls `singlerot` with `overtone_a` and `qnmr`, then plots `real(spectrum)` with `plot_1d`.

Although the header mentions an unphysically strong pulse, this file assigns no explicit `rf_pwr`, `rf_dur`, `rf_frq`, or `method` field. It therefore documents no numerical pulse setting; the executable setup shown is the `singlerot`/`overtone_a` call and the panoramic sweep. This avoids transferring the separate overtone_pa average-treatment RF settings from mas_boron_2.m to this example.
