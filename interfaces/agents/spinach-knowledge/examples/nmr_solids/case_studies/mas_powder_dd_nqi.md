# examples/nmr_solids/case_studies/mas_powder_dd_nqi.m

- Signature: `mas_powder_dd_nqi()`

## Purpose

Powder magic angle spinning spectrum of a pair of dipole-coupled quadrupolar nuclei; this is apparently something that other simu- lation packages cannot do. Parameters from Jeongjae Lee. Calculation time: seconds

## Physical / mathematical content
- Simulates a powder magic-angle-spinning spectrum of dipole-coupled quadrupolar ²³Na and ¹⁷O nuclei at 9.4 T, using coordinates and electric-field-gradient-derived quadrupolar interactions.
- Acquires the ¹⁷O signal at 100 kHz spinning over a 200-point spherical orientation grid; applies exponential apodisation and a zero-filled Fourier transform before plotting.

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
