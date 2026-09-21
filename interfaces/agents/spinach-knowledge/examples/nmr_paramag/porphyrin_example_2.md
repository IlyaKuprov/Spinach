# examples/nmr_paramag/porphyrin_example_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/porphyrin_example_2.m`
- Signature: `porphyrin_example_2()`
- Total lines: 76

## Purpose

Computing PCS using different models in a basic Cu(II) porphyrin complex. The metal is at the origin. See the "getting started" manual at The paper describing the distributed PCS model used below is available at

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

## Implementation structure

- Computing PCS using different models in a basic Cu(II) porphyrin complex.
- The metal is at the origin. See the "getting started" manual at
- The paper describing the distributed PCS model used below is available at
- Porphyrin ring proton coordinates
- Cu(II) g-tensor
- Curie susceptibility tensor
- Metal position
- PCS calculation using the point model
- Coordinates of all atoms with significant spin population
- Mulliken spin populations
- Multipole ranks to include
- Multipole moments

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2chi()`, `ppcs()`, `points2mult()`, `lpcs()`, `oparse()`, `dft_pcs_cu()`, `hfc2pcs()`.
