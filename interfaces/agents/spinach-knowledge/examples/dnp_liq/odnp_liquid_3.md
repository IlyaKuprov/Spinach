# examples/dnp_liq/odnp_liquid_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/odnp_liquid_3.m`
- Signature: `odnp_liquid_3()`
- Total lines: 95

## Purpose

Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment with an electron and a nucleus connected by a hyperfine coupling. A g-hyperfine cross-correla- tion effect is visible at high field. The steady state is computed by setting the time derivative to zero in the inhomogeneous master equation, and solving the resulting algebraic equation for the steady s

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Steady state nuclear magnetisation as a function of microwave frequency
- offset and the magnet field in a DNP experiment with an electron and a
- nucleus connected by a hyperfine coupling. A g-hyperfine cross-correla-
- tion effect is visible at high field.
- The steady state is computed by setting the time derivative to zero in
- the inhomogeneous master equation, and solving the resulting algebraic
- equation for the steady state density matrix.
- Calculation time: minutes.
- Spin system
- Anisotropic Zeeman interactions
- Isotropic hyperfine coupling
- Coordinates for dipolar coupling

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `field_grid()`, `create()`, `basis()`, `state()`, `operator()`, `answer()`, `liquid()`, `kfigure()`, `xlim()`, `kxlabel()`, `kylabel()`, `kcolourbar()`.
