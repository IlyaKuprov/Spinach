# examples/nmr_paramag/porphyrin_example_3.m

- Signature: `porphyrin_example_3()`

## Purpose

Computes PCS using different models in basic Cu(II) porphyrin complex. See the "getting started" manual at The paper describing the distributed PCS model used below is available at Calculation time: minutes, 64GB of RAM required

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

## Implementation structure

- Computes PCS using different models in basic Cu(II) porphyrin complex. See
- the "getting started" manual at
- The paper describing the distributed PCS model used below is available at
- Calculation time: minutes, 64GB of RAM required
- Porphyrin ring proton coordinates
- Cu(II) g-tensor eigenvalues
- Curie susceptibility tensor
- Parse ORCA log
- Extract hyperfine tensors
- Compute HFC PCS
- Parse ORCA cube and pad the density with zeros to avoid PBC effects
- Compute extentens for original density without padding
