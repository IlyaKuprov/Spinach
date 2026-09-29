# experiments/pseudocon/pcs_combi_fit.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/pcs_combi_fit.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=pcs_combi_fit.m)

## Purpose

Searches ambiguous diamagnetic and paramagnetic shift assignments, scoring each pairing by how well a fitted susceptibility tensor explains the resulting pseudocontact shifts. The routine performs assignment enumeration and PCS fitting; it is not a pulse-sequence simulation.

## Inputs

The input is a `parameters` structure with these fields:

- `hfcs`: cell array of 3-by-3 hyperfine tensors in Gauss, typically from `gparse`.
- `isotopes`: cell array of isotope-name character strings corresponding to the hyperfine tensors, for example {'1H','1H'}.
- `spin_groups`: cell array of integer index vectors, one per distinct chemical-shift entry, listing the hyperfine tensors/spins associated with that shift; the source gives `{[28 22 30]; [25 33 81]}` as an example.
- `d_shifts` and `p_shifts`: vectors of distinct diamagnetic and paramagnetic chemical shifts, respectively, in ppm.
- `d_ambig` and `p_ambig`: cell arrays of integer index groups marking shift entries whose diamagnetic or paramagnetic assignment may be permuted. These indices rearrange the corresponding shift vectors.

The source checks required fields and rejects an obsolete `nel` field because hyperfine tensors are normalised per unpaired electron. Expanded hyperfine tensors, isotope labels, and shift arrays are passed to `pcs2chi`, which checks their compatible lengths and tensor/isotope types.

## Assignment search and outputs

For each ambiguity group, `perms` enumerates the possible orderings; the routine forms the Cartesian product of these permutations independently for the diamagnetic and paramagnetic shifts. For every pair of complete assignments, it subtracts the diamagnetic shift vector from the paramagnetic vector. The per-group PCS values are expanded over `spin_groups`, pairing each spin's `hfcs` and isotope, and passed to `pcs2chi`. The assignment score is the `pcs2chi` least-squares error; the minimum-scoring pair is selected. Candidate scoring uses `parfor`.

The selected per-group shifts are returned as `d_shifts` and `p_shifts`. The routine then returns per-spin predicted PCS values from `hfc2pcs`, selected per-spin experimental PCS values (paramagnetic minus diamagnetic), and fitted susceptibility `chi`. `total_theo` adds predicted PCS values to the selected diamagnetic shifts repeated according to the number of spins in each `spin_groups` entry.

The underlying PCS tensor fit uses the model cited by `pcs2chi`: [DOI 10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).