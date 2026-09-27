# examples/relaxation_theory/trosy_fluorine_sym.m

- Signature: `trosy_fluorine_sym()`

## Purpose

Calculate transverse relaxation and broad and narrow TROSY line rates as functions of magnetic field for `19F` and its directly bonded `13C` in a 3-fluorotyrosine-labelled protein. The calculation also separates the `13C` TROSY rate into dipole–dipole (DD), chemical-shift anisotropy (CSA), and DD–CSA cross-correlation contributions. The source describes the analytical calculation time as seconds.

## Physical / mathematical content

- The calculation uses the fluorine and carbon shielding tensors and atomic coordinates extracted from a 3-fluorotyrosine DFT calculation.
- `rlx_dd_csa` returns total `R2` rates and broad and narrow TROSY rates for both nuclei, plus the DD, CSA, and cross-correlation contributions used for the `13C` mechanism plot.

## Numerical / algorithmic content

- Read `../standard_systems/3_fluoro_tyr.log` with `gparse` and `g2spinach`, mapping carbon to `13C` and fluorine to `19F` with the supplied values `[186.38 192.97]`.
- Extract shielding tensors and coordinates from entries 8 (`19F`) and 7 (`13C`) of the resulting interaction data.
- Evaluate 20 proton Larmor frequencies from 200 to 800 MHz and convert them to magnetic fields using `spin('1H')`.
- At each field, call `rlx_dd_csa` with a `25e-9` s timescale, the nuclei `{'19F','13C'}`, their shielding tensors, and their coordinates.

## Implementation structure

- Plot broad TROSY, total `R2`, and narrow TROSY rates against proton Larmor frequency separately for `19F` and `13C`.
- Plot the `13C` TROSY DD and CSA contributions alongside `-abs` of the DD–CSA cross-correlation contribution in a stacked bar chart. Plot rates are labelled in Hz.