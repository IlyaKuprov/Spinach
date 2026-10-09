# examples/nmr_overtone/mas_boron_1.m

- MATLAB implementation: [examples/nmr_overtone/mas_boron_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/mas_boron_1.m)

This example simulates a 10B magic-angle-spinning (MAS) overtone spectrum in Z-detection. The source comments say the sample spins in the JEOL direction, credit parameters to Nghia Duong and Yusuke Nishiyama, and identify the target as the most intense of five overtone spinning sidebands. The comment estimates hours of calculation. These comments do not assert an experimental fit or reproduction.

The system has isotope 10B and magnet setting 16.4. Its only listed coupling is `eeqq2nqi(0.7e6,0.0,3,[0 0 0])`; the source does not state units for the magnet or coupling arguments. The relaxation settings are `damp`, diagonal retention, zero equilibrium, and `damp_rate=50`. The basis is `sphten-liouv` with approximation `none`.

The sequence uses rank 12, axis `[sqrt(2/3) 0 sqrt(1/3)]`, rate 70000, and grid `rep_2ang_800pts_sph`. It sets sweep `[-141e3 -139e3]`, 256 points, 256-point zero-fill, and `axis_units='kHz'`. Both initial state and receiver are 10B `Lz`, matching the source's Z-detection description. The source does not set explicit RF power, duration, frequency, or an average-treatment field in this file.

The simulation call is `singlerot` with `overtone_a` and `qnmr`. The script plots `real(spectrum)` through `plot_1d`; it does not apply a separate phase factor. Its narrow sweep and spherical 800-point grid distinguish this sideband-focused setup from the two panoramic 10B examples in this group.
