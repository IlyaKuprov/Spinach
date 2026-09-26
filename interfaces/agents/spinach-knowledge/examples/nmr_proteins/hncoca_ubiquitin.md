# examples/nmr_proteins/hncoca_ubiquitin.m

- Signature: `hncoca_ubiquitin()`

## Purpose

Theoretical HN(CO)CA of human ubiquitin. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Theoretical HN(CO)CA of human ubiquitin. It is assumed that
- only the backbone is 13C,15N-labelled.
- Calculation time: minutes, faster with a Tesla A100 GPU.
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
