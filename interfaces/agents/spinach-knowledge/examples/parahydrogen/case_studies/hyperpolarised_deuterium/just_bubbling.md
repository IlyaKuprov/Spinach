# examples/parahydrogen/case_studies/hyperpolarised_deuterium/just_bubbling.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/case_studies/hyperpolarised_deuterium/just_bubbling.m`
- Signature: `just_bubbling()`
- Total lines: 115

## Purpose

Evolution of state populations under ortho-deuterium bubbling in the presence of a parahydrogenation cata- lyst. No pulses, just evolution. Paper link to follow in due course.

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Evolution of state populations under ortho-deuterium
- bubbling in the presence of a parahydrogenation cata-
- lyst. No pulses, just evolution. Paper link to follow
- in due course.
- Spin system
- Experimental chemical shifts
- Experimental J-couplings
- NQI tensors from a DFT calculation
- Cartesian coordinates, DFT calculation
- Kinetics
- Magnet field
- Simulation formalsim

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxan()`, `deut_pair()`, `unit_state()`, `coils()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `pumped_state()`, `magpump()`, `evolution()`, `traj_a()`, `kfigure()`, `scale_figure()`.
