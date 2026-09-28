# examples/nmr_paramag/porphyrin_example_2.m

- Signature: `porphyrin_example_2()`

## Purpose

Compares point-model, distributed-density, and DFT-hyperfine PCS for a Cu(II) porphyrin complex, with the metal at the origin. The source links to the [pseudocontact-shift analysis getting-started manual](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis) and the [paper describing the distributed PCS model](http://dx.doi.org/10.1039/c6cp05437d).

## Physical / mathematical content

The example compares PCS predicted by a point electron, a distributed electron density represented by atomic coordinates and Mulliken spin populations, and DFT hyperfine tensors. The distributed calculation expands the density in multipoles.

## Numerical / algorithmic content

The Cu(II) tensor is derived from `g_cu=diag([2.0000 2.0000 2.2000])` at 298 K for spin 1/2. The distributed model uses Mulliken populations `[0.6 0.1 0.1 0.1 0.1]` and multipole ranks `0:14`; the DFT comparison uses hyperfine tensors 26:37 from `cu_porph_hfc.out`.

## Implementation structure

- Define porphyrin proton coordinates and compute point-model PCS at the origin.
- Build the spin-population multipole representation and calculate distributed PCS.
- Parse the ORCA output, calculate HFC-derived PCS, and display the three-way comparison.
