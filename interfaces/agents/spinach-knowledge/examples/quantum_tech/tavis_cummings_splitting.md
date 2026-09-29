# examples/quantum_tech/tavis_cummings_splitting.m

Source: [examples/quantum_tech/tavis_cummings_splitting.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/tavis_cummings_splitting.m)

- Signature: `tavis_cummings_splitting()`

## Model

This is a closed, resonant Tavis–Cummings calculation for one to four identical electron spins coupled to one cavity mode. The isotope list contains `E` for each electron and `C3` for the cavity truncation; it does not specify a diamond defect or a measured sample. The magnetic field and both rotating-frame mode frequencies are zero. Every electron is coupled to the shared cavity with `coupling=3e6`, or 3 MHz in the frequency convention used for the plot. Thus the model has zero detuning, no applied drive, and no configured dissipation.

The full Zeeman-Hilbert basis is built without approximation, then the cavity-context Hamiltonian is restricted to the one-excitation sector. The selected basis states use electron Zeeman labels `ZL1` and `ZL2` with the cavity states `BL1` and `BL2`: one state has one electron flipped and the cavity unexcited, and the remaining state has all spins unflipped and one cavity excitation. This is the spin-flip/cavity-photon selection used here; there is no defect-specific EPR spectrum, g tensor, hyperfine structure, field sweep, or powder/orientation average.

## Observable and plot

For each spin count, the script diagonalises the projected Hamiltonian and takes the span between its lowest and highest eigenvalues as the bright-mode splitting. It compares that gap with `2*sqrt(n)*2*pi*coupling` and raises an error if the relative discrepancy exceeds `1e-10`. The plot places spin count on the horizontal axis and splitting in MHz on the vertical axis, overlaying numerical points and the square-root prediction. It depicts the ideal collective coupling model, not an experimental spectrum.

The source cites Tavis and Cummings, *Physical Review* **170**, 379 (1968) ([DOI](https://doi.org/10.1103/PhysRev.170.379)).
