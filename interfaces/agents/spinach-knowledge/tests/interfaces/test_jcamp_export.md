# tests/interfaces/test_jcamp_export.m

`result=test_jcamp_export()` exercises the JCAMP exporter on synthetic data
with explicit coordinates and amplitudes. It checks descending and irregular
axes, floating-point precision, complex channels, ragged multidimensional
pages, differing page grids, assigned peaks, EMR peak descriptors, linked
NMR/EPR blocks, missing observations, record wrapping, byte-identical file output, and refusal of
invalid inputs (including malformed delays, isotope labels, units, and EMR
multiplicity, ordinate units, methods, and labels) without overwriting an
existing destination. Repeated page coordinates remain distinct, and JCAMP
5.01 magnitude/power unit labels do not recompute caller data and are rejected
for complex NMR traces. NMR assignments cover optional heights, widths,
multiplicity, and method comments.

This regression test is not a spectrometer-vendor acceptance test. Its
uncompressed pair reader checks the samples emitted by the exporter; external
JCAMP readers are needed to assess reader-specific compatibility.

- [MATLAB test](https://github.com/IlyaKuprov/Spinach/blob/main/tests/interfaces/test_jcamp_export.m)
- [Exporter](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jcamp/jcamp_export.m)
