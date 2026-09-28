# examples/visualisation/cst_peptide_bond.m

- Signature: `cst_peptide_bond()`

## Purpose

Read the alanine Gaussian log at `../standard_systems/amino_acids/ala.log` and visualise the carbon, proton, and nitrogen shielding tensors. The source notes that antisymmetric shielding-tensor components are ignored.

## Implementation

The script parses the log with `gparse`, then makes a three-panel figure using `cst_display` in harmonic style. It selects C, H, and N in turn, passing display parameters `0.01`, `0.05`, and `0.01`, respectively.
