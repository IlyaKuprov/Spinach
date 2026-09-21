# examples/giant_spin/case_studies/mn3_trimer_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_levels.m`
- Signature: `mn3_trimer_levels()`
- Total lines: 63

## Purpose

Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, from zero to 10 Tesla in the full 216-state Hilbert space. The lowest levels are the ones that the 16-state and 26-state effective bases of the paper are built to reproduce. Reproduces Figure 5 of https://arxiv.org/abs/2609.16352. Calculation time: seconds.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Zeeman energy level diagram of the (CH6N3)2MnCl4 molecular crystal, a
- linear trimer of three S=5/2 manganese ions with isotropic exchange
- between neighbours and an axial plus rhombic zero-field splitting on
- every ion, from zero to 10 Tesla in the full 216-state Hilbert space.
- The lowest levels are the ones that the 16-state and 26-state effec-
- tive bases of the paper are built to reproduce. Reproduces Figure 5
- Calculation time: seconds
- Magnet must be 1 Tesla, the field is set below
- Parallel pool size
- Three S=5/2 spins with g=2
- Isotropic exchange, J=-2.42 cm^-1 in the H=-2*J*S1*S2 convention of the paper
- Zero-field splitting, D=0.167 cm^-1 and E=0.040 cm^-1 on every ion

## Internal Spinach / MATLAB structure cues

- Called routines in the main body: `icm2hz()`, `zfs2mat()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `orientation()`, `eig()`, `sort()`, `hz2icm()`, `kfigure()`, `plot()`, `kxlabel()`, `kylabel()`, `save()`.
