# etc/gui/IO Formats/IOSpinSystem.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/IOSpinSystem.cpp)

The common reader accepts an expected field count, column indices and a coordinate scale. Those are format-specific: a Cartesian table in bohr needs a different scale from an XYZ row in angstroms. Parsing into an existing frame updates the atom records rather than concatenating successive geometries.

The matrix reader extracts three specified columns from each of three rows, normalizes Fortran D exponents and checks field counts before indexing. A malformed first numeric column is not a zero tensor. Zero quadrupolar tensors are skipped, but a zero final row alone does not make the whole matrix zero.
