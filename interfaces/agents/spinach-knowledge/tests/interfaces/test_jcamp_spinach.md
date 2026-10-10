# tests/interfaces/test_jcamp_spinach.m

Tests Spinach-facing JCAMP wrappers across 1D pulse acquisition, odd/even FFT
axes, scalar homonuclear settings, 2D States and cosine/sine components,
3D four-branch FIDs and processed tensors, mixed domains, field/ENDOR scans,
HYSCORE grids, and explicit pulse delays. Numerical table readback checks
coordinates and signed amplitudes; metadata checks cover observation frequency,
isotope order, page coordinates, and named component identities. Ambiguous
shapes and non-nuclear NMR observations must be refused.

[Source](../../../../../tests/interfaces/test_jcamp_spinach.m)
