# interfaces/jcamp/jcamp_signal.m

Exports explicitly sampled Spinach EMR results with
`text=jcamp_signal(spin_system,axes,signal,units,names,info)`.
One to three physical axis columns are supplied in MATLAB array order,
with a unit and name per axis. Signal dimensions must match the coordinates;
a 1D result is a column. Named component structures are also accepted.
This accommodates irregular DEER/ESEEM delays, ENDOR scans, and custom
multidimensional maps without an inferred acquisition grid.

Ownership, filename, detection mode, EMR core method, physical parameter
description, and extra metadata are explicit in info. The model field comes
from spin_system in tesla (the fieldsweep scaling field is not a measured field), simulation source is Spinach, and data type is EMR SIMULATION.
Units use the JCAMP EMR vocabulary, e.g. SECOND, HERTZ, and TESLA. The
internal file structure and pages are built by jcamp_grid and written by
jcamp_export. No sorting, interpolation, signal processing, or phase cycling
occurs. Standard NMR outputs use jcamp_nmr instead.

[Source](../../../../jcamp/jcamp_signal.m) · [Usage examples](../../../../jcamp/README.md#export-directly-from-spinach-results)
