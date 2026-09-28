# examples/kinetics/glucose_exsy_a.m

- Signature: `glucose_exsy_a()`

## Purpose

2D EXSY of transmembrane exchange of 2,2,3,3-tetrafluoroglucose. See the fitting example set for the script that yielded the parameters used below. Calculation time: seconds.

## Physical / mathematical content

Sixteen ¹⁹F spins represent the α and β forms on the inside and outside of the membrane (four four-spin subsystems). Their chemical shifts, couplings, and coordinates are specified separately. A four-state reaction-rate matrix describes translocation; the initial population vector passed to `equilibrate` has an α/β imbalance. Relaxation is Redfield with secular retention and separate correlation times for inside and outside states.

## Numerical / algorithmic content

The script runs a 2D `noesy` EXSY simulation with a 0.5 s mixing time, then applies squared-cosine apodisation and States-style quadrature processing before the two Fourier transforms. It loads the experimental spectrum, ranks it, and plots the simulated and experimental spectra plus a pointwise-deviation histogram.

## Implementation structure

- Uses B₀ = 9.4 T and ¹⁹F chemical shifts/J couplings for the four α/β, inside/outside subsystems.
- Sets the kinetic starting vector to [3.2258, 0, 3.1902, 0] before equilibration.
- Uses `parameters.sweep=[8000 8000]`, `npoints=[256 256]`, and `zerofill=[1024 512]` (ppm axes).
- Loads `glucose_expt_a.mat` / `Expression1` for the experimental comparison.
