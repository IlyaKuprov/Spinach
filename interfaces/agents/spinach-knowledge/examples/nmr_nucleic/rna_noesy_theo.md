# examples/nmr_nucleic/rna_noesy_theo.m

- Signature: `rna_noesy_theo()`

## Purpose

1H-1H NOESY spectrum of the example RNA molecule provided by the Gerhard Wagner group at Harvard University. Calculation time: hours Shunsuke Imai Scott Robson Gerhard Wagner Zenawi Welderufael Ilya Kuprov

## Physical / mathematical content

- Nucleic-acid NMR examples. These files specialise biomolecular NMR workflows to RNA or DNA systems, with labelled nuclei, residue-level assignments, and multidimensional heteronuclear transfer logic.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- 1H-1H NOESY spectrum of the example RNA molecule provided by
- the Gerhard Wagner group at Harvard University.
- Calculation time: hours
- Shunsuke Imai
- Scott Robson
- Gerhard Wagner
- Zenawi Welderufael
- Ilya Kuprov
- Import RNA data
- Magnet field
- Tolerances
- Relaxation theory
