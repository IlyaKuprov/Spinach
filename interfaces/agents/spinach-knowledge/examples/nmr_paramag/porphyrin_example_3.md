# examples/nmr_paramag/porphyrin_example_3.m

- Signature: `porphyrin_example_3()`
- Source: [examples/nmr_paramag/porphyrin_example_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/porphyrin_example_3.m)

## Purpose

Compares proton pseudocontact shifts (PCS) for a Cu(II) porphyrin using hyperfine couplings and a distributed spin-density calculation. The source cites the [PCS analysis manual](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis) and the [distributed PCS paper](http://dx.doi.org/10.1039/c6cp05437d). The source estimates minutes of calculation time and 64 GB of RAM.

## Physical / mathematical content

- Defines 12 porphyrin-ring proton coordinates and a diagonal Cu(II) g tensor with entries 2.0000, 2.0000, and 2.2000. It converts this tensor to a Curie susceptibility with `g2chi(g_cu,298,1/2)` (298 K, spin 1/2).
- Reads `cu_porph_hfc.out` with `oparse`, takes hyperfine tensors 26 through 37, and calculates the proton PCS from each tensor with `hfc2pcs`.
- Reads the spin-density cube `cu_porph_sd120.spindens.3d` with `ocparse`, pads it with zeros using `pad_size=2`, then calls `kpcs` in `fft` mode to calculate the density-based PCS.

## Numerical / algorithmic content

The displayed comparison is the HFC and PDE PCS values in ppm. The script also plots the spin-density and PCS volumes with the parsed molecular geometry. It reads ORCA calculation outputs; it does not load a measured NMR spectrum or establish agreement with an experiment. The source does not annotate tensor or coordinate units, so none are assigned here.

## Implementation structure

- Defines the proton coordinates and Cu(II) g tensor, then forms the susceptibility at 298 K for spin 1/2.
- Parses the ORCA hyperfine output and spin-density cube, and compares discrete HFC-derived PCS with PCS from the padded density through the Kuprov equation.
- Displays the numerical shifts and two molecular-volume schematics.
