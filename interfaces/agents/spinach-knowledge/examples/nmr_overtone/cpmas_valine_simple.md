# examples/nmr_overtone/cpmas_valine_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/cpmas_valine_simple.m`
- Signature: `cpmas_valine_simple()`
- Total lines: 75

## Purpose

Cross-polarization experiment between protons and 14N overtone transition in N-acetylvaline under MAS, computed using Fokker- Planck formalism. Valine quadrupolar tensor data comes from: Calculation time: hours.

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Cross-polarization experiment between protons and 14N overtone
- transition in N-acetylvaline under MAS, computed using Fokker-
- Planck formalism. Valine quadrupolar tensor data comes from:
- Calculation time: hours.
- System specification
- Relaxation theory
- Basis set
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Spectrum setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `singlerot()`, `kfigure()`, `plot_1d()`.
