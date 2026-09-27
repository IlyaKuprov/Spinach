# examples/nmr_proteins/noesyhsqc_ubiquitin_prot.m

- Signature: `noesyhsqc_ubiquitin_prot()`

## Purpose

1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with 65 ms mixing time. It is assumed that the protein is not 13C-labelled. Calculation time: hours.

## Physical / mathematical content
- Simulates a 3D ¹H–¹H–¹⁵N NOESY-HSQC spectrum of ¹⁵N-labelled, non-¹³C-labelled ubiquitin from `1D3Z.pdb` and `1D3Z.bmrb` at 900 MHz, with a 65 ms mixing time.
- Uses Redfield relaxation with a 5 ns correlation time and a 90 Hz coupling parameter; the simulated dimensions have 128 × 64 × 128 points.
- Applies squared-cosine apodisation to four signal components, then Fourier transforms the three dimensions with zero filling to 512 × 256 × 512 points.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900
- MHz with 65 ms mixing time. It is assumed that the protein is
- not 13C-labelled.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure
- Kill carbons
