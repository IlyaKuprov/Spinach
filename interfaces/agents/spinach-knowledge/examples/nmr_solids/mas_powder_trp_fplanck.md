# examples/nmr_solids/mas_powder_trp_fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_trp_fplanck.m`
- Signature: `mas_powder_trp_fplanck()`
- Total lines: 106

## Purpose

13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Isotropic chemical shifts come from the experimental data. Coordinates are from X-ray data and CSAs are estimated with DFT. Calculation time: days, hours with a Tesla A100 GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H),
- computed using the Fokker-Planck MAS formalism. Isotropic chemical
- shifts come from the experimental data. Coordinates are from X-ray
- data and CSAs are estimated with DFT.
- Calculation time: days, hours with a Tesla A100 GPU.
- % First molecule in the unit cell
- Spin system properties (DFT calculation)
- Magnet field
- Experimental chemical shifts, first conformation
- Basis set
- Algorithmic options
- sys.enable={'gpu'};

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
