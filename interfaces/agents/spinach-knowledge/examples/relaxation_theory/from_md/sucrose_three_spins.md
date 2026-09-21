# examples/relaxation_theory/from_md/sucrose_three_spins.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/from_md/sucrose_three_spins.m`
- Signature: `sucrose_three_spins()`
- Total lines: 112

## Purpose

One of the calculations reported in the JMR paper with Jim Prestegard: a three-spin subsystem from the glucose ring of sucrose -an illustra- tion of incorrect viscosity of TIP3P water. The experimental value of tau_c for sucrose is around 90 ps (and this is correctly reproduced by OPC and TIP5P water), but TIP3P only agrees with Redfield theory when tau_c is set to 37 ps in the latter. Here, numerical relaxation supe

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- One of the calculations reported in the JMR paper with Jim Prestegard:
- a three-spin subsystem from the glucose ring of sucrose -an illustra-
- tion of incorrect viscosity of TIP3P water.
- The experimental value of tau_c for sucrose is around 90 ps (and this
- is correctly reproduced by OPC and TIP5P water), but TIP3P only agrees
- with Redfield theory when tau_c is set to 37 ps in the latter.
- Here, numerical relaxation superoperator computed from a long MD tra-
- jectory is compared with the analytical one computed using the isotro-
- pic rotational diffusion approximation.
- Calculation time: minutes, with most of the time spent
- computing MD frame Hamiltonians
- Three-spin system

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `load()`, `traj()`, `report()`, `double()`, `traj_slice()`, `dipolar()`, `orientation()`, `ngce()`, `relaxation()`, `kfigure()`, `errorbar()`, `kxlabel()`.
