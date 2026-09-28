# examples/nmr_paramag/porphyrin_example_1.m

- Signature: `porphyrin_example_1()`

## Purpose

Computes point-model PCS for porphyrin-ring protons in basic Cu(II) and Co(II) complexes, with the metal at the origin. The source links to the [pseudocontact-shift analysis getting-started manual](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Physical / mathematical content

The calculation derives Curie susceptibility tensors from the Cu(II) and Co(II) g-tensors, then evaluates the point-model pseudocontact shifts at the listed proton coordinates.

## Numerical / algorithmic content

The two tensors are calculated at 298 K for spin 1/2: `g_co=diag([3.0 3.0 2.0])` and `g_cu=diag([2.0 2.0 2.2])`. The metal coordinate is `[0 0 0]`; `ppcs` returns separate PCS values for Co and Cu.

## Implementation structure

- Define the porphyrin proton coordinates and the two metal-ion g-tensors.
- Convert the g-tensors to Curie susceptibility tensors at 298 K for spin 1/2.
- Calculate and display the Co(II) and Cu(II) point-model PCS values.
