# Pulsed-field magnetometry case studies

Spinach reproductions of the simulation figures of the qdmag paper (Liu, Chen, Cupo, Fry, Cheng, arXiv:2609.16352), which describes non-equilibrium magnetisation of molecular magnets under a time-varying field with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. Every script uses the paper's own Hamiltonian parameters, temperatures, spectral density, sweep profiles, and stair widths; the kernel functions involved are `stevens.m` and `stev2sph.m` (crystal field to rank 12 through `inter.giant`), `phonon_oper.m` (the thermally dressed coupling operator, rebuilt on every stair and applied as Hilbert space matrix products; `rlx_phonon.m` is the standalone Liouville space form of the same dissipator and is not called here), and `pulsed_field.m` (the experiment, run under `crystal` or `powder`).

| Script | Paper figure | Content |
| --- | --- | --- |
| `ho_pzdo4_profiles.m` | Figure 2 | Ho(pzdo)4, J=8, four field profiles: linear, piecewise linear, spline of a measured 65 T pulse, sinusoidal at the clock gap |
| `ho_pzdo4_powder.m` | Figure 3 | Ho(pzdo)4, single-orientation curves on a Lebedev grid, powder average, and equilibrium |
| `dimer_exchange_types.m` | Figure 4 | Two S=1/2 spins with isotropic, anisotropic, and antisymmetric exchange, against equilibrium |
| `mn3_trimer_levels.m` | Figure 5 | (CH6N3)2MnCl4 trimer of three S=5/2, Zeeman level diagram in the full 216-state space |
| `mn3_trimer_magn.m` | Figure 6 | The same trimer, 50 T/ms sweep at 0.6 K in the full space, against equilibrium |

`ho_pzdo4_params.m` holds the 90 crystal-field coefficients shared by the two Ho scripts. The trimer scripts use the full Hilbert space where the paper uses 16-state and 26-state effective bases; at 216 states no reduction is needed. Each script saves its curves to a `.mat` file next to it.
