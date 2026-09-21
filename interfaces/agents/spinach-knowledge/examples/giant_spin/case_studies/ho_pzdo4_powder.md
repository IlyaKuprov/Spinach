# examples/giant_spin/case_studies/ho_pzdo4_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_powder.m`
- Signature: `ho_pzdo4_powder()`
- Total lines: 98

## Purpose

Powder-averaged pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework, a J=8 giant spin with a crystal field to twelfth spherical rank, under a 10 T/ms linear sweep at 2 K with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The magnetisation of every orientation of the two-angle Lebedev grid is plotted alongside the powder average and the thermal equilibrium magnetisation. Reproduces Fig 3 of https://arxiv.org/abs/2609.16352. Calculation time: hours.

## Physical / mathematical content

- One `E17` spin (J=8) with the effective g-factor 1.24 and the 90 CASSCF crystal field coefficients of `ho_pzdo4_params.m`, ranks 2 to 12, converted rank by rank from Stevens to spherical tensor form with `icm2hz` and `stev2sph` into `inter.giant.coeff`.
- Every orientation of the grid stays below the equilibrium curve over the whole sweep, and so does the powder average, as the paper states: at 2 K the phonon bath cannot follow a 10 T/ms sweep.
- The spin-phonon coupling operator has unit elements between adjacent m_J states; super-Ohmic bath (`phonon_alpha=2`) with lambda=10 cm^-1 and I0=1e-14 ps/rad; the observable is the moment -1.24*J_z in Bohr magnetons.

## Numerical / algorithmic content

- `zeeman-hilb`, `powder` context with the `labframe` assumption, `parameters.needs={'zeeman_op'}`, `parameters.sum_up=false` so that the answer of every orientation is returned, and the `leb_2ang_rank_29` grid; `sys.magnet=1`; 10 ns stairs, 10^5 stairs, one record every 1000 stairs; 16 workers. J_z comes from `operator(spin_system,'Lz','E17')`.
- The powder average is the weight-sum of the single-orientation observables. The equilibrium curve is the same weight-sum of `real(hdot(coil,equilibrium(spin_system,H)))`, with H the Hilbert space Hamiltonian of each grid orientation (isotropic part plus `orientation` of the anisotropic part, minus the unit-field Zeeman operator, plus the field times the Zeeman operator per Tesla) at every recorded field.
