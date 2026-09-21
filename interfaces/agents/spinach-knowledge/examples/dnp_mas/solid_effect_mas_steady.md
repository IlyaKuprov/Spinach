# examples/dnp_mas/solid_effect_mas_steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/solid_effect_mas_steady.m`
- Signature: `solid_effect_mas_steady()`
- Total lines: 108

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Steady state rotor period simulation for a single crystal, computed using Newton-Raphson steady state solver. Calculation time: seconds

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A MAS DNP simulation performed as described in Fred Mentink-
- Vigier's paper (Spinach rotation conventions are different):
- Steady state rotor period simulation for a single crystal,
- computed using Newton-Raphson steady state solver.
- Calculation time: seconds
- Magnet field
- Spin specification
- Interactions
- Relaxation parameters
- Basis set
- Spinach housekeeping
- Stack generation parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rotor_stack()`, `operator()`, `relaxation()`, `speye()`, `propagator()`, `steady()`, `rho()`, `step()`, `kfigure()`, `trajan()`, `equilibrium()`, `hamiltonian()`, `assume()`, `state()`.
