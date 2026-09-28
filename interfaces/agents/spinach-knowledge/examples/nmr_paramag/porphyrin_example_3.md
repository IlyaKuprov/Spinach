# examples/nmr_paramag/porphyrin_example_3.m

- Signature: `porphyrin_example_3()`

## Purpose

Computes pseudocontact shifts (PCS) using hyperfine tensors and a distributed spin-density model for a Cu(II) porphyrin complex. See the [PCS analysis manual](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis) and the [paper describing the distributed PCS model](http://dx.doi.org/10.1039/c6cp05437d). Calculation time: minutes, 64GB of RAM required.

## Physical / mathematical content

- Converts the Cu(II) g-tensor to a Curie susceptibility tensor at 298 K for spin 1/2.
- Computes proton PCS from hyperfine tensors and solves the Kuprov equation using the padded spin density and susceptibility tensor.

## Numerical / algorithmic content

- Uses FFT mode for the density-based PCS calculation and displays the HFC and PDE PCS values side by side in ppm.

## Implementation structure

- Defines porphyrin ring proton coordinates and Cu(II) g-tensor eigenvalues.
- Parses an ORCA log for hyperfine tensors and an ORCA spin-density cube with zero padding to avoid PBC effects.
- Plots the spin density and PCS field alongside the molecular geometry.
