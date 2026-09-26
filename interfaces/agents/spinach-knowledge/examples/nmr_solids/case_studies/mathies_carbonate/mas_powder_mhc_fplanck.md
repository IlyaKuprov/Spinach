# examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck.m

- Signature: `mas_powder_mhc_fplanck()`

## Purpose

All protons in the unit cell of monohydrocalcite, magic angle spinning NMR simulation. Further details in: Calculation time: hours, minutes on a GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- All protons in the unit cell of monohydrocalcite, magic
- angle spinning NMR simulation. Further details in:
- Calculation time: hours, minutes on a GPU.
- 400 MHz NMR
- Read CASTEP file
- Drop C, O, and Ca atoms
- Get isotopes
- Convert shielding tensors into shift using the
- parametrisation of Huang et al. ACIE 2021
- Get coordinates
- Basis set
- Interaction cut-off, Hz
