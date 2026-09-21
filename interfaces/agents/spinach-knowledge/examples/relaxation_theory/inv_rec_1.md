# examples/relaxation_theory/inv_rec_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/inv_rec_1.m`
- Signature: `inv_rec_1()`
- Total lines: 63

## Purpose

A simple inversion-recovery experiment; longitudinal magnetisation is monitored as a function of time. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A simple inversion-recovery experiment; longitudinal
- magnetisation is monitored as a function of time.
- Calculation time: seconds.
- Spin system
- Zeeman interactions
- Complete basis set
- Relaxation theory
- Spinach housekeeping
- Isotropic thermal equilibrium
- Detection state
- Static Liouvillian superoperator
- Pulse operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `state()`, `hamiltonian()`, `assume()`, `operator()`, `relaxation()`, `step()`, `evolution()`, `kfigure()`, `scale_figure()`, `kylabel()`, `kxlabel()`.
