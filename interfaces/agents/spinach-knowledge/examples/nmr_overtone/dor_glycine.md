# examples/nmr_overtone/dor_glycine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/dor_glycine.m`
- Signature: `dor_glycine()`
- Total lines: 78

## Purpose

Panoramic double rotation overtone 14N spectrum of glycine, simulated as described in our paper (Figure 1B): A short pulse with instrumentally inaccessible power is gi- ven to make the excitation pattern uniform. Glycine quadru- polar tensor data comes from O'Dell and Ratcliffe: Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Panoramic double rotation overtone 14N spectrum of glycine,
- simulated as described in our paper (Figure 1B):
- A short pulse with instrumentally inaccessible power is gi-
- ven to make the excitation pattern uniform. Glycine quadru-
- polar tensor data comes from O'Dell and Ratcliffe:
- Calculation time: hours
- System specification
- Relaxation theory
- Basis set
- Algorithmic options
- Spinach housekeeping
- Magic angle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `doublerot()`, `kfigure()`, `plot_1d()`.
