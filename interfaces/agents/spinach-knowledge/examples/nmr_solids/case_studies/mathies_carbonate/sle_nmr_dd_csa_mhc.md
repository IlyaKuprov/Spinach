# examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m`
- Signature: `sle_nmr_dd_csa_mhc()`
- Total lines: 101

## Purpose

Water protons in the unit cell of monohydrocalcite, inc- luding slow isotropic rotational diffusion and MAS. Fur- ther details in: Calculation time: minutes, seconds with a GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Water protons in the unit cell of monohydrocalcite, inc-
- luding slow isotropic rotational diffusion and MAS. Fur-
- ther details in:
- Calculation time: minutes, seconds with a GPU.
- 400 MHz NMR
- Read CASTEP file
- Drop C, O, and Ca atoms
- Keep two protons
- Convert shielding tensors into shift using the
- parametrisation of Huang et al. ACIE 2021
- Get coordinates
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `ismember()`, `create()`, `basis()`, `state()`, `kfigure()`, `tau_c()`, `max_rank()`, `gridfree()`, `apodisation()`, `fftshift()`, `plot_1d()`, `ylim()`, `klegend()`, `kylabel()`.
