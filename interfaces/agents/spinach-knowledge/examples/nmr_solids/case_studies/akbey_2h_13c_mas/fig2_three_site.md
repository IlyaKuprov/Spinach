# examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m

- Signature: `fig2_three_site()`

## Purpose

Three-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the quadru- polar tensor. Set to reproduce Figure 2 in: Calculation time: seconds.

## Physical / mathematical content
- Simulates three-site deuterium position exchange at 9.4 T; the sites differ in chemical shift and quadrupolar-tensor orientation, have equal populations, and exchange pairwise at 10⁴ s⁻¹.
- Acquires a magic-angle-spinning spectrum at 8.5 kHz using an 800-point spherical orientation grid, then Fourier-transforms the 1,024-point signal without apodisation.

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
