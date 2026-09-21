# examples/parahydrogen/altadena_propanal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/altadena_propanal.m`
- Signature: `altadena_propanal()`
- Total lines: 68

## Purpose

ALTADENA experiment simulation for the parahydrogenation of acrolein into propanal. Simple model of the ALTADENA effect is used: perfect- ly adiabatic transfer is assumed and the isotropic mixing in low fi- eld is ignored completely. Note the small flip angle. Calculation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- ALTADENA experiment simulation for the parahydrogenation of acrolein
- into propanal. Simple model of the ALTADENA effect is used: perfect-
- ly adiabatic transfer is assumed and the isotropic mixing in low fi-
- eld is ignored completely. Note the small flip angle.
- Calculation time: seconds
- Spin system
- Magnetic field
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
