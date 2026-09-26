# examples/nmr_overtone/dante_glycine.m

- Signature: `dante_glycine()`

## Purpose

14N overtone DANTE spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tensor data comes from the paper by O'Dell and Ratcliffe: Calculation time: minutes

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- 14N overtone DANTE spectrum of glycine, computed using Fokker-Planck
- formalism. Glycine quadrupolar tensor data comes from the paper by
- O'Dell and Ratcliffe:
- Calculation time: minutes
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Spectrum setup
- Simulation
