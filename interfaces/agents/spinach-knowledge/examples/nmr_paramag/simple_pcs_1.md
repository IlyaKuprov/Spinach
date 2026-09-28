# examples/nmr_paramag/simple_pcs_1.m

- Signature: `simple_pcs_1()`

## Purpose

Pseudocontact shift and Curie relaxation on a proton due to the presence of a point magnetic susceptibility centre. Calculation time: seconds

## Physical / mathematical content

- Defines a susceptibility tensor and its position relative to a proton at the origin.
- Specifies Redfield relaxation with zero equilibrium, lab-frame terms, and a 10 ps correlation time.

## Numerical / algorithmic content

- Creates a spin system in the `sphten-liouv` formalism without basis approximation, then calculates R1 and R2 from the relaxation superoperator and proton `Lz` and `L+` states.

## Implementation structure

- Sets the 14.1 T magnetic field, proton isotope, diamagnetic shift, coordinates, susceptibility, and relaxation parameters.
- Builds the spin system and basis, calculates relaxation rates, and displays R1 and R2.
