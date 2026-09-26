# examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m

- Signature: `fig2_three_site()`

## Purpose

Three-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the quadru- polar tensor. Set to reproduce Figure 2 in: Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Three-site position exchange for a deuterium nucleus. The sites
- differ in the chemical shift and the orientation of the quadru-
- polar tensor. Set to reproduce Figure 2 in:
- Calculation time: seconds.
- Magnet field
- Spin system
- Quadrupolar interactions
- Kinetics
- Basis set
- Spinach housekeeping
- Sequence parameters
- MAS parameters
