# examples/singlet_states/warren_singlet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/warren_singlet.m`
- Signature: `warren_singlet()`
- Total lines: 45

## Purpose

A demonstration that long-lived states exist that are immune not only to dipolar and CSA, but also to quadrupolar relaxati- on in certain circumstances. Full Redfield superoperator for dipolar and quadrupolar relaxation in liquid state is compu- ted and diagonalized. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- A demonstration that long-lived states exist that are immune
- not only to dipolar and CSA, but also to quadrupolar relaxati-
- on in certain circumstances. Full Redfield superoperator for
- dipolar and quadrupolar relaxation in liquid state is compu-
- ted and diagonalized.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Relaxation superoperator accuracy
- Basis set
- Spinach housekeeping
- Relaxation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `relaxation()`.
