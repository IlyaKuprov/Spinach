# examples/relaxation_theory/trosy_fluorine_sym.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/trosy_fluorine_sym.m`
- Signature: `trosy_fluorine_sym()`
- Total lines: 74

## Purpose

Transverse relaxation rate as a function of the applied magnetic field in a 3-fluorotyrosine labelled protein. The fluorine atom and its directly bonded carbon are included. Analytical calcula- tions broken down by mechanism. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

## Implementation structure

- Transverse relaxation rate as a function of the applied magnetic
- field in a 3-fluorotyrosine labelled protein. The fluorine atom
- and its directly bonded carbon are included. Analytical calcula-
- tions broken down by mechanism.
- Calculation time: seconds.
- Read 3-fluorotyrosine DFT calculation
- Extract coordinates and CSAs
- Magnetic field grid
- Loop over magnetic fields
- Call the analytical function
- Relaxation rates
- Mechanisms for 13C

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `spin()`, `rlx_dd_csa()`, `r2c()`, `r2f()`, `f_bro()`, `f_nar()`, `c_bro()`, `c_nar()`, `c_tro_dd()`, `c_tro_csa()`, `c_tro_xc()`, `kfigure()`, `ylim()`, `kxlabel()`.
