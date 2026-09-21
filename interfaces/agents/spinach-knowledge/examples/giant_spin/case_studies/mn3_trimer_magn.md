# examples/giant_spin/case_studies/mn3_trimer_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_magn.m`
- Signature: `mn3_trimer_magn()`
- Total lines: 88

## Purpose

Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The full 216-state Hilbert space is used; the paper solves the same problem in 

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a
- linear trimer of three S=5/2 manganese ions with isotropic exchange
- between neighbours and an axial plus rhombic zero-field splitting on
- every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon
- relaxation in the generalised Lindblad form of Saito and Miyashita.
- The full 216-state Hilbert space is used; the paper solves the same
- problem in 16-state and 26-state effective bases. The thermal equi-
- librium magnetisation is plotted for comparison. Reproduces Figure 6
- Calculation time: minutes
- Magnet must be 1 Tesla, the field is set by the sweep
- Parallel pool size
- Three S=5/2 spins with g=2, phonon bath at 0.6 K

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `zfs2mat()`, `create()`, `basis()`, `operator()`, `double()`, `crystal()`, `hamiltonian()`, `assume()`, `orientation()`, `m_eq()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`, `save()`.
