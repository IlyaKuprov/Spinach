# examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m`
- Signature: `mas_powder_gly_13c()`
- Total lines: 103

## Purpose

13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism and a spherical grid. The field dependence of the line shape of 13CA due to the presence of the quadrupolar 14N nucleus is shown. The calcula- tion is performed in the rotating frame with respect to 13C and the laboratory frame with respect to 14N. Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C MAS spectrum of glycine powder (assuming decoupling of 1H),
- computed using the Fokker-Planck MAS formalism and a spherical
- grid. The field dependence of the line shape of 13CA due to the
- presence of the quadrupolar 14N nucleus is shown. The calcula-
- tion is performed in the rotating frame with respect to 13C and
- the laboratory frame with respect to 14N.
- Calculation time: seconds.
- Read CASTEP file
- Drop H and O atoms
- Keep 13CA and 14N
- Convert shielding tensors into shift
- Set isotropic chemical shifts to experimental values

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `shift_iso()`, `castep2nqi()`, `remtrace()`, `kfigure()`, `fields()`, `create()`, `basis()`, `spin()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `plot_1d()`, `klegend()`.
