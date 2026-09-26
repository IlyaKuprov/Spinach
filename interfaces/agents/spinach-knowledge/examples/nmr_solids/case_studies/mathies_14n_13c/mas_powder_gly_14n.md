# examples/nmr_solids/case_studies/mathies_14n_13c/powder_gly_14n.m

- Signature: `powder_gly_14n()`

## Purpose

Static 14N powder spectrum of glycine (assuming decoupling of 1H and 13C), computed on a spherical grid from CASTEP tensors and compared to a simulation using the quadrupolar parameters measured by O'Dell and Schurko (PCCP 2009, DOI 10.1039/b906114b). Numerical rotating frame transformation is used because 14N quadrupolar interaction is large. Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static conditions: chemical-shift anisotropy, quadrupolar coupling, and orientation averaging using direct powder quadrature on a spherical grid.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Static 14N powder spectrum of glycine (assuming decoupling of
- 1H and 13C), computed on a spherical grid from CASTEP tensors
- and compared to a simulation using the quadrupolar parameters
- measured by O'Dell and Schurko (PCCP 2009, DOI 10.1039/b906114b).
- Numerical rotating frame transformation is used because 14N
- quadrupolar interaction is large.
- Calculation time: seconds.
- Read CASTEP file
- Drop H, O, and C atoms
- keep only 1 14N
- Convert shielding tensor into shift
- Set isotropic chemical shift to experimental value
- Quadrupolar interaction from CASTEP
- Magnet field
