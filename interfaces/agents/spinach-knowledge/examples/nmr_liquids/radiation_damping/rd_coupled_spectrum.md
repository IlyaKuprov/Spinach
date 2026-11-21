# rd_coupled_spectrum: radiation damping

[MATLAB example](../../../../../../examples/nmr_liquids/radiation_damping/rd_coupled_spectrum.m)

Simulates an inverted, strongly coupled proton pair with and without radiation
damping. The feedback is determined by the total proton transverse
magnetisation; its phase and strength are editable example parameters.
Cartesian magnetisation traces accompany magnitude spectra of the two FIDs.

Spectral processing uses an exponential window, the first-point correction,
and zero filling. The no-feedback calculation is the internal comparison;
this script is a fixed-parameter demonstration, not an interactive interface.
