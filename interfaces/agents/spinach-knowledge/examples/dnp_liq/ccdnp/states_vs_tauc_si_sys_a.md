# examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_a.m`
- Signature: `states_vs_tauc_si_sys_a()`
- Total lines: 143

## Purpose

Steady state populations of various spin states a function of rotational correlation time in a DNP experiment with two electrons connected by ex- change coupling, both coupled to a nucleus by dipolar couplings. Further particulars in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Steady state populations of various spin states a function of rotational
- correlation time in a DNP experiment with two electrons connected by ex-
- change coupling, both coupled to a nucleus by dipolar couplings. Further
- particulars in:
- Calculation time: seconds
- Magnetic field, Tesla
- Spin system
- Zeeman interactions
- Exchange coupling
- Coordinates for anisotropic HF
- Basis set
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `tau_c()`, `create()`, `basis()`, `state()`, `operator()`, `answer()`, `liquid()`, `kfigure()`, `subplot()`, `scale_figure()`, `kxlabel()`, `kylabel()`, `xlim()`.
