# interfaces/jcamp/jcamp_nmr.m

Exports native Spinach 1D, 2D, and 3D NMR arrays or scalar structures of named
quadrature components. Syntax: `text=jcamp_nmr(spin_system,parameters,signal,domains,info)`.
The `domains` cell row declares time/frequency in physical F1/F2/F3 order.
Arrays match example layouts: 1D columns, 2D [F2,F1], 3D [F1,F2,F3].
Sequence sweep widths and offsets are Hz; scalar sweep/offset/nucleus entries
are shared across physical dimensions. Actual array lengths determine sample
counts. Time axes start at acquisition zero with dwell 1/sweep. Frequency axes
use `ft_axis`, requiring at least three points per frequency dimension.

The field and isotopes determine observation frequencies in MHz; axes stay
in Hz regardless of plot display units. Named components become separate
LINK blocks; complex arrays retain signed real and imaginary pages. The
writer input is built internally and passed to `jcamp_export`. No FFT,
conjugation, phase cycling, or quadrature recombination occurs. Mixed-domain
blocks take the type of their tabulated first-array-dimension axis.

`info` supplies ownership, filename (empty for text only), pulse-sequence text,
and additional metadata. Time-containing arrays also require the actual
JCAMP delay pair in microseconds and acquisition convention. The wrappers
record all axis nuclei, observation frequencies, offsets, and domains in
private metadata. Use `jcamp_grid` plus the low-level writer for non-standard
NMR sampling; use `jcamp_epr` for electron results.

[Source](../../../../jcamp/jcamp_nmr.m) · [Usage examples](../../../../jcamp/README.md#export-directly-from-spinach-results)
