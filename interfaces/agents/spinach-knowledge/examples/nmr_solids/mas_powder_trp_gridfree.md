# examples/nmr_solids/mas_powder_trp_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_trp_gridfree.m`
- Signature: `mas_powder_trp_gridfree()`
- Total lines: 118

## Purpose

13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism. Isotro- pic chemical shifts come from the experimental data. Coordinates and CSAs are estimated with DFT. A polyadic representation of the evolution generator is used, further particulars here: Calculation time: hours on a Tesla V100 GPU, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H),
- computed using the grid-free Fokker-Planck MAS formalism. Isotro-
- pic chemical shifts come from the experimental data. Coordinates
- and CSAs are estimated with DFT. A polyadic representation of the
- evolution generator is used, further particulars here:
- Calculation time: hours on a Tesla V100 GPU,
- much longer on CPU
- % First molecule in the unit cell
- Spin system properties (DFT calculation)
- Magnet field
- First conformation
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `gridfree()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
