# examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m

- Signature: `mas_powder_gly_13c()`

## Purpose

13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism and a spherical grid. The field dependence of the line shape of 13CA due to the presence of the quadrupolar 14N nucleus is shown. The calcula- tion is performed in the rotating frame with respect to 13C and the laboratory frame with respect to 14N. Calculation time: seconds.

## Physical / mathematical content
- Simulates glycine’s ¹³C magic-angle-spinning powder spectrum, assuming ¹H decoupling, with the Fokker–Planck MAS formalism and a spherical orientation grid.
- Retains ¹³Cα and quadrupolar ¹⁴N from CASTEP data to examine the ¹⁴N interaction’s effect on the ¹³Cα line shape at 4.7, 9.4, and 14.1 T.
- Acquires the ¹³C signal at 10 kHz spinning, then applies exponential apodisation and a zero-filled Fourier transform to plot the spectra.

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
