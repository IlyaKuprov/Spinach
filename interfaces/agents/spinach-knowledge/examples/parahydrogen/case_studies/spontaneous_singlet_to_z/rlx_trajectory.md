# examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m`
- Signature: `rlx_trajectory()`
- Total lines: 75

## Purpose

Time dependence of LzSz and Lz+Sz spin orders in a para- hydrogen molecule coordinated to a nickel cage that cre- ates large chemical shift anisotropy. Done for Gloggler group, paper link coming in due course; see also the Ma- thematica worksheet.

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Time dependence of LzSz and Lz+Sz spin orders in a para-
- hydrogen molecule coordinated to a nickel cage that cre-
- ates large chemical shift anisotropy. Done for Gloggler
- group, paper link coming in due course; see also the Ma-
- thematica worksheet.
- Magnet field, Tesla
- Isotopes
- Coordinates (for dipole tensor)
- Zeeman interaction tensors, traces subtracted
- Relaxation theory parameters
- Formalism and basis
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `relaxation()`, `singlet()`, `state()`, `evolution()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`, `scale_figure()`.
