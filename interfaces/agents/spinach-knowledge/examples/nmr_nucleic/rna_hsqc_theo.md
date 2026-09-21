# examples/nmr_nucleic/rna_hsqc_theo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_nucleic/rna_hsqc_theo.m`
- Signature: `rna_hsqc_theo()`
- Total lines: 79

## Purpose

1H-13C HSQC spectrum of the example RNA molecule provided by the Wagner group. Calculation time: minutes Shunsuke Imai Scott Robson Gerhard Wagner Zenawi Welderufael Ilya Kuprov

## Physical / mathematical content

- Nucleic-acid NMR examples. These files specialise biomolecular NMR workflows to RNA or DNA systems, with labelled nuclei, residue-level assignments, and multidimensional heteronuclear transfer logic.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- 1H-13C HSQC spectrum of the example RNA molecule provided by
- the Wagner group.
- Calculation time: minutes
- Shunsuke Imai
- Scott Robson
- Gerhard Wagner
- Zenawi Welderufael
- Ilya Kuprov
- Import RNA data
- Magnet field
- Tolerances
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `nuclacid()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
