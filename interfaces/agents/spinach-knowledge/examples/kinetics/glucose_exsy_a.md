# examples/kinetics/glucose_exsy_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/glucose_exsy_a.m`
- Signature: `glucose_exsy_a()`
- Total lines: 168

## Purpose

2D EXSY of transmembrane exchange of 2,2,3,3-tetrafluoroglucose. See the fitting example set for the script that yielded the parameters used below. Calculation time: seconds

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 2D EXSY of transmembrane exchange of 2,2,3,3-tetrafluoroglucose. See
- the fitting example set for the script that yielded the parameters
- used below.
- Calculation time: seconds
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Cartesian coordinates
- Chemical subsystems
- Reaction rate matrix
- Equilibrate translocation with alpha-beta imbalance as the start

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `equilibrate()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `load()`, `rot90()`, `keep_rank()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_2d()`, `ktitle()`, `histogram()`.
