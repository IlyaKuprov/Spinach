# examples/nmr_overtone/powder_glycine_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/powder_glycine_2.m`
- Signature: `powder_glycine_2()`
- Total lines: 60

## Purpose

Overtone detection 14N powder NMR spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tensor data comes from the paper by O'Dell and Ratcliffe: This simulation demonstrates that the spin state that gives rise to the overtone signal in a static sample is the T2,-2 coherence. Calculation time: seconds

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Overtone detection 14N powder NMR spectrum of glycine, computed using
- Fokker-Planck formalism. Glycine quadrupolar tensor data comes from
- the paper by O'Dell and Ratcliffe:
- This simulation demonstrates that the spin state that gives rise to
- the overtone signal in a static sample is the T2,-2 coherence.
- Calculation time: seconds
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `powder()`, `kfigure()`, `plot_1d()`.
