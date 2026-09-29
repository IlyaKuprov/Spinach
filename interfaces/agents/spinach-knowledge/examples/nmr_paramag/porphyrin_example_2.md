# examples/nmr_paramag/porphyrin_example_2.m

- Signature: `porphyrin_example_2()`
- Source: [examples/nmr_paramag/porphyrin_example_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/porphyrin_example_2.m)
- Manual cited by the source: [Pseudocontact shift analysis](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis)
- Distributed-model paper cited by the source: [DOI 10.1039/c6cp05437d](https://doi.org/10.1039/c6cp05437d)

## Cu(II) porphyrin comparison

Twelve hard-coded porphyrin-ring proton coordinates are evaluated relative to the metal at `[0 0 0]`. The Cu(II) g-tensor is `diag([2.0000 2.0000 2.2000])`; `g2chi(g_cu,298,1/2)` supplies the Curie susceptibility tensor. The source does not label the units of the coordinate values or the `298` argument.

The point result is calculated with `ppcs`. For the distributed model, the five spin-bearing coordinates (as listed in the source) are `[0,0,0]`, `[-1.452280156504,1.452280156496,0]`, `[1.452280156504,-1.452280156496,0]`, `[-1.452280156496,-1.452280156504,0]`, and `[1.452280156496,1.452280156504,0]`. Their Mulliken spin populations are `[0.6,0.1,0.1,0.1,0.1]`. The source does not state coordinate units. It obtains moments for ranks 0 through 14 with `points2mult` and calculates PCS using `lpcs`.

A third column is derived from precomputed hyperfine tensors: `oparse('cu_porph_hfc.out')` reads the ORCA log, entries 26 through 37 are selected, and `hfc2pcs(...,chi_cu,'1H')` converts them to proton PCS. The script prints point, distributed, and HFC-derived PCS columns in ppm; it does not run the ORCA calculation.

## Scope and omissions

This is a basic Cu(II) porphyrin example, not a carbonic-anhydrase case. It has no protein residue/site, measured spectrum, spectral simulation, specified magnetic field, or field/temperature sweep. The source computes three PCS comparisons from supplied coordinates and a supplied ORCA output file; it performs no parameter-fitting call.
