# examples/nmr_solids/case_studies/mas_powder_dd_nqi.m

- Signature: `mas_powder_dd_nqi()`

## Purpose

Powder magic angle spinning spectrum of a pair of dipole-coupled quadrupolar nuclei; this is apparently something that other simu- lation packages cannot do. Parameters from Jeongjae Lee. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder magic angle spinning spectrum of a pair of dipole-coupled
- quadrupolar nuclei; this is apparently something that other simu-
- lation packages cannot do. Parameters from Jeongjae Lee.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Enable GPU
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup
- Simulation
