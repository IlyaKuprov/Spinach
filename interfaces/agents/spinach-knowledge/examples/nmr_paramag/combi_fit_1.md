# examples/nmr_paramag/combi_fit_1.m

- Signature: `combi_fit_1()`

## Purpose

Extracting the susceptibility tensor from DFT hyperfine tensors and experimental paramagnetic shifts. A combinatorial procedure is used that cycles through ambiguous assignments. Calculation time: minutes.

## Physical / mathematical content

- The routine enumerates the specified ambiguous diamagnetic- and paramagnetic-shift assignments and fits a susceptibility tensor from the DFT hyperfine tensors to the experimental paramagnetic shifts.

## Numerical / algorithmic content

The setup supplies 27 proton isotope labels, nine spin groups, measured diamagnetic and paramagnetic shifts, and the corresponding ambiguity sets to `pcs_combi_fit`. The returned theoretical and experimental PCS values are plotted against one another.

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
