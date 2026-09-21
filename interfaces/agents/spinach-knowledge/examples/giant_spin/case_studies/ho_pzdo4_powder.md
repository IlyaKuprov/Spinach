# examples/giant_spin/case_studies/ho_pzdo4_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_powder.m`
- Signature: `ho_pzdo4_powder()`
- Total lines: 98

## Purpose

Powder-averaged pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework, a J=8 giant spin with a crystal field to twelfth spherical rank, under a 10 T/ms linear sweep at 2 K with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The magnetisation of every orientation of the two-angle Lebedev grid is plotted alongside the powder average and the thermal equilibrium magnetisation. Reproduces Fig 3 of https://arxiv.org/abs/2609.16352. Calculation time: hours.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The powder-averaged thermal equilibrium magnetisation at every recorded field is obtained from `equilibrium()` applied to the field-dependent Hilbert space Hamiltonian of each grid orientation, read out with `hdot()`, and weighted with the grid weights.

## Implementation structure

- Powder-averaged pulsed-field magnetisation of the Ho(pzdo)4 metal-
- organic framework, a J=8 giant spin with a crystal field to twelfth
- spherical rank, under a 10 T/ms linear sweep at 2 K with spin-phonon
- relaxation in the generalised Lindblad form of Saito and Miyashita.
- The magnetisation of every orientation of the two-angle Lebedev grid
- is plotted alongside the powder average and the thermal equilibrium
- magnetisation. Reproduces Fig 3 of
- Calculation time: hours
- Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention
- Convert Stevens coefficients into spherical tensor coefficients, Hz, rank by rank
- Magnet must be 1 Tesla, the field is set by the sweep
- Parallel pool size

## Internal Spinach / MATLAB structure cues

- Called routines in the main body: `ho_pzdo4_params()`, `stev2sph()`, `icm2hz()`, `create()`, `basis()`, `operator()`, `double()`, `powder()`, `hamiltonian()`, `assume()`, `orientation()`, `equilibrium()`, `hdot()`, `kfigure()`, `plot()`.
