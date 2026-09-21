# examples/nmr_paramag/combi_fit_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/combi_fit_1.m`
- Signature: `combi_fit_1()`
- Total lines: 50

## Purpose

Extracting the susceptibility tensor from DFT hyperfine tensors and experimental paramagnetic shifts. A combinatorial procedure is used that cycles through ambiguous assignments. Calculation time: minutes.

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

## Implementation structure

- Extracting the susceptibility tensor from DFT hyperfine tensors and
- experimental paramagnetic shifts. A combinatorial procedure is used
- that cycles through ambiguous assignments.
- Calculation time: minutes.
- Read DFT HFCs in Gauss
- Isotope list
- Spin groups with identical PCS
- Diamagnetic shifts
- Diamagnetic shift ambiguities
- Paramagnetic shifts
- Paramagnetic shift ambiguities
- Run the combinatorial fitting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `pcs_combi_fit()`, `kfigure()`.
