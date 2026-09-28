# examples/kinetics/glucose_exsy_b.m

- Signature: `glucose_exsy_b()`

## Purpose

2D EXSY of transmembrane exchange of 3,3-difluoroglucose. See the fitting example set for the script that yielded the parameters used below. Calculation time: seconds.

## Physical / mathematical content

Eight ¹⁹F spins form four two-spin α/β, inside/outside subsystems, with separate shifts, couplings, and coordinates. A four-state reaction-rate matrix models translocation, and the source equilibrates a starting population vector with an α/β imbalance. Relaxation combines Redfield and T1/T2 terms with secular retention and state-specific correlation times.

## Numerical / algorithmic content

The 2D `noesy` simulation uses a 0.5 s mixing period, squared-cosine apodisation, States processing, and Fourier transforms in both dimensions. The processed calculation is compared with the loaded experimental spectrum; the script also plots their pointwise deviation histogram.

## Implementation structure

- Uses B₀ = 9.3933 T and the four two-spin α/β, inside/outside ¹⁹F subsystems.
- Sets the starting concentration vector to [3.8034, 0, 14.2442, 0] for equilibration.
- Uses `parameters.sweep=[8650 8650]`, `npoints=[512 1024]`, and `zerofill=[1024 1024]` (ppm axes).
- Loads `glucose_expt_b.mat` / `spec` for the experimental comparison.
