# kernel/utilities/isot2elem.m

`elements=isot2elem(isotopes)` strips digits from each character row vector in
a cell array, preserving its shape and order. For example,
`isot2elem({'1H','13C','35Cl'})` returns `{'H','C','Cl'}`. Empty cell arrays
retain their dimensions, and unsupported input types are rejected.

This is a literal string conversion, not a physical-isotope lookup: all
non-digit characters remain unchanged. In particular, an isomer suffix in
`'99Tc_m'` remains in the output `'Tc_m'`. Use ordinary mass-number-plus-element
labels when element symbols are required.
