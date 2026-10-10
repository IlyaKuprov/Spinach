# interfaces/jcamp/jcamp_grid.m

`blocks=jcamp_grid(axes,signal,units,names,block)` is the shared sampled-array
to JCAMP-block translator. It does not write a file. Axis columns, units,
and names are in MATLAB dimension order; one to three dimensions are
supported. A 1D signal is a column. Higher-dimensional array sizes must
match every axis length. A scalar structure of arrays represents named
quadrature/receiver components, retained as separate blocks rather than
recombined. The supplied block contains title, type, and technique metadata.

Real 1D signals remain traces; complex traces are split by the final writer.
For multidimensional signals, dimension 1 is tabulated and remaining physical
coordinates label pages. Real and imaginary ordinates are separately named,
with original signs and tensor ordering intact. All component arrays must
have the same coordinate shape. The helper supports explicit irregular NMR
sampling when paired with appropriate NMR metadata and jcamp_export; the
public EMR explicit-axis entry point is jcamp_signal.

[Source](../../../../jcamp/jcamp_grid.m) · [Usage examples](../../../../jcamp/README.md#export-directly-from-spinach-results)
