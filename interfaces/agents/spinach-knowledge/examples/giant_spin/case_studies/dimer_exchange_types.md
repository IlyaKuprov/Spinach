# examples/giant_spin/case_studies/dimer_exchange_types.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/dimer_exchange_types.m`
- Signature: `dimer_exchange_types()`
- Total lines: 96

## Purpose

Pulsed-field magnetisation of a dimer of two S=1/2 spins with four types of exchange coupling tensor: isotropic, two anisotropic, and antisymmetric, at 0.2 K under a 10 T/ms sweep to 1 T, with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The out-of-equilibrium curves are compared with the thermal equilibrium magnetisation. Reproduces Figure 4 of https://arxiv.org/abs/2609.16352. Calculation time: minutes.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Pulsed-field magnetisation of a dimer of two S=1/2 spins with four
- types of exchange coupling tensor: isotropic, two anisotropic, and
- antisymmetric, at 0.2 K under a 10 T/ms sweep to 1 T, with spin-
- phonon relaxation in the generalised Lindblad form of Saito and
- Miyashita. The out-of-equilibrium curves are compared with the
- thermal equilibrium magnetisation. Reproduces Figure 4 of
- Calculation time: minutes
- Magnet must be 1 Tesla, the field is set by the sweep
- Parallel pool size
- Two electron spins with g=2
- Temperature of the phonon bath
- Formalism and basis set

## Internal Spinach / MATLAB structure cues

- Called routines in the main body: `stevens()`, `double()`, `kfigure()`, `scale_figure()`, `icm2hz()`, `create()`, `basis()`, `crystal()`, `hamiltonian()`, `assume()`, `orientation()`, `eig()`, `trace()`, `subplot()`, `plot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
