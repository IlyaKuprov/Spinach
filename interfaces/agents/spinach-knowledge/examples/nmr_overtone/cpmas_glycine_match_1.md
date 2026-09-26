# examples/nmr_overtone/cpmas_glycine_match_1.m

- Signature: `cpmas_glycine_match_1()`

## Purpose

Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tensor da- ta comes from the paper by O'Dell and Ratcliffe: Hartmann-Hahn condition profile with a rough powder grid, as a function of 1H RF power. Calculation time: minutes

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Cross-polarization experiment between protons and 14N overtone
- transition in glycine under MAS. Glycine quadrupolar tensor da-
- ta comes from the paper by O'Dell and Ratcliffe:
- Hartmann-Hahn condition profile with a rough powder grid, as a
- function of 1H RF power.
- Calculation time: minutes
- System specification
- Relaxation theory
- Basis set
- Spinach housekeeping
- Magic angle
- Spectrum setup
