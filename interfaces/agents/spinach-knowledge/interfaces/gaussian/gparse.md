# interfaces/gaussian/gparse.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/gaussian/gparse.m) · [Spinach Wiki: gparse.m](https://spindynamics.org/wiki/index.php?title=gparse.m)

- Signature: `props=gparse(filename,options)`

## Input and options

`filename` must be a nonempty character vector naming an existing Gaussian 03, 09, or 16 output log. The optional `options` argument is a cell array of character vectors; valid entries are `'g_nosymm'`, `'cst_nosymm'`, and `'hfc_nosymm'`. By default, parsed g, chemical-shielding, and hyperfine tensors are symmetrised; each corresponding option disables that symmetrisation. Other entries and non-cell/non-character options are rejected.

The parser reads the existing text log using MATLAB file and text-processing functions; it does not launch Gaussian. The source recommends a route section containing `#p nmr=(giao,spinspin,susceptibility) output=pickett pop=minimal IOp(6/82=1)` to produce a useful property log. This is guidance, not an argument checked by the parser.

## Returned property structure

The single output is a structure, not a Spinach system. Parsed fields include input and standard geometries in Å, atom count, method, SCF energy in Hartree, atomic numbers and symbols, charge, multiplicity, spin components and `s_sq`, and isotope mass numbers, atomic masses, nuclear spins, quadrupole moments, and magnetic moments when present. Special geometry symbols include `Bq` for ghost atoms, `X` for dummy atoms, and `TV` for translation vectors. If no standard orientation is printed but an input geometry exists, `std_geom` is set to that input geometry.

Magnetic-property fields are conditional on the printed sections: `hfc.iso` and `hfc.full.matrix` are hyperfine values/tensors in Gauss; `g_tensor.matrix` and its eigenvalues/eigenvectors describe the dimensionless g tensor; `cst` contains absolute shielding tensors; `k_couplings` and `j_couplings` are isotropic couplings in Hz; `srt` contains spin-rotation tensors in Hz; and `nqi` contains nuclear quadrupolar tensors in Hz. Other parsed values include susceptibility, Gibbs free energy in Hartree, and electric dipole moment in Debye. `filename`, `error`, and `complete` identify the source log and termination state.

## Transformations and checks

For Link1/multi-job logs, the parser processes the occurrences in the text and retains the latest occurrence of each parsed property rather than returning one structure per job. It reads geometries and the Gaussian isotope table, and constructs atomic symbols from the atomic numbers. The anisotropic spin-dipole part of the HFC is normalised by `multiplicity-1` before it is combined with the isotropic Fermi-contact term; closed-shell HFC sections are ignored. If an anisotropic HFC section appears without multiplicity, parsing stops with an error.

The parser checks the HFC unit columns using the factor `2.802495`, checks g-tensor eigenvalues against printed g shifts referenced to `g_e=2.0023193043`, and checks shielding trace/3 against the printed isotropic value (tolerance `1e-3`). Disagreements raise errors rather than silently returning the inconsistent tensor. Gaussian error-termination markers set `props.error`; normal termination sets `props.complete`. An incomplete log without a detected error produces a warning, so inspect these fields before downstream use.

[Spinach Wiki: gparse.m](https://spindynamics.org/wiki/index.php?title=gparse.m)
