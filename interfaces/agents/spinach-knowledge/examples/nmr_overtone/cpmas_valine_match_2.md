# examples/nmr_overtone/cpmas_valine_match_2.m

- Signature: `cpmas_valine_match_2()`

## Purpose

Cross-polarization experiment between protons and 14N overtone transition in N-acetylvaline under MAS, computed using Fokker- Planck formalism. Valine quadrupolar tensor data comes from: Hartmann-Hahn condition profile with a rough powder grid, as a function of 1H RF power and spinning rate. Calculation time: hours.

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Cross-polarization experiment between protons and 14N overtone
- transition in N-acetylvaline under MAS, computed using Fokker-
- Planck formalism. Valine quadrupolar tensor data comes from:
- Hartmann-Hahn condition profile with a rough powder grid, as a
- function of 1H RF power and spinning rate.
- Calculation time: hours.
- System specification
- Relaxation theory
- Basis set
- Spinach housekeeping
- Magic angle
- Spectrum setup
