# experiments/pseudocon/pcs_combi_fit.m

- Signature: `[d_shifts,p_shifts,pcs_theo,pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)`

## Purpose

Combinatorial PCS fitting function. Takes into account potential am- biguities in diamagnetic and paramagnetic NMR assignments. Syntax: [d_shifts,p_shifts,pcs_theo,... pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- parameters.hfcs -cell array of 3x3 hyperfine tensors,
- in Gauss, usually out of gparse() or
- something similar
- parameters.isotopes -cell array of isotope specificati-
- ons, e.g. {'1H','1H'}
- parameters.spin_groups -cell array of integer vectors
- specifying the numbers of spins
- that have each of the chemical
- shifts specified, e.g.
- {[28 22 30]; [25 33 81]}
- parameters.d_shifts -a vector of unique diamagnetic che-
- mical shifts, in ppm
- parameters.p_shifts -a vector of unique paramagnetic
- chemical shifts, in ppm
- parameters.d_ambig -a cell array of integer vectors spe-
- cifying the spins for which the dia-
- magnetic assignment can potentially
- be swapped around.
- parameters.p_ambig -a cell array of integer vectors spe-
- cifying the spins for which the para-
- magnetic assignment can potentially
- be swapped around.

## Outputs

- d_shifts -diamagnetic chemical shifts, optimally permuted
- p_shifts -paramagnetic chemical shifts, optimally permuted
- pcs_theo -theoretical pseudocontact shifts
- pcs_expt -experimental pseudocontact shifts from optimally
- permuted assignments
- chi -rank 2 part of the magnetic susceptibility tensor,
- in cubic Angstrom
- total_theo -theoretical total NMR chemcial shifts, computed
- as a sum of d_shifts and pcs_theo

## Implementation structure

- Combinatorial PCS fitting function. Takes into account potential am-
- biguities in diamagnetic and paramagnetic NMR assignments. Syntax:
- [d_shifts,p_shifts,pcs_theo,...
- pcs_expt,chi,total_theo]=pcs_combi_fit(parameters)
- parameters.hfcs -cell array of 3x3 hyperfine tensors,
- in Gauss, usually out of gparse() or
- something similar
- parameters.isotopes -cell array of isotope specificati-
- ons, e.g. {'1H','1H'}
- parameters.spin_groups -cell array of integer vectors
- specifying the numbers of spins
- that have each of the chemical
