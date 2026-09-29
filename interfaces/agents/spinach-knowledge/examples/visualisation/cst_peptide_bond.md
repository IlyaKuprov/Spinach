# examples/visualisation/cst_peptide_bond.m

- Signature: `cst_peptide_bond()`
- Source: [examples/visualisation/cst_peptide_bond.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/visualisation/cst_peptide_bond.m)

## Purpose and input

This example parses the Gaussian log `../standard_systems/amino_acids/ala.log` with `gparse` and visualises shielding tensors for the alanine peptide-bond example. The source comment states that antisymmetric components of the shielding tensors are ignored.

## Rendering

A three-panel figure uses spherical-harmonic rendering throughout. The calls select C, H, and N and pass display parameters 0.01, 0.05, and 0.01, respectively; the panels are titled ¹³C CST, ¹H CST, and ¹⁵N CST. Each panel sets camera position [40 40 40], and the figure uses `scale_figure([2.0 1.0])`. The source gives no units or further interpretation for the display parameters, so they are recorded as the values supplied to `cst_display` rather than physical tensor magnitudes.
