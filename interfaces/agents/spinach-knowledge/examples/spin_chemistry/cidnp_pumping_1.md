# examples/spin_chemistry/cidnp_pumping_1.m

- Signature: `cidnp_pumping_1()`

## Purpose

A simulation of the matrix in Equation 2 of IK's paper on chemically amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011). Calculation time: seconds.

## Physical / mathematical content


## Physical / mathematical content

- At 14.1 T, the model is a 1H-19F pair with zero isotropic shifts, 50 Hz scalar coupling, and a fluorine CSA tensor with principal values `[-47 -16 63]`. The spins are separated by `2.60` along y.
- The Redfield relaxation model retains secular terms, uses `tau_c=110e-12`, and is configured for equilibrium magnetisation at `298 K`. The script adds pumping terms with strengths `1.3` for 1H and `34.0` for 19F.

## Numerical / algorithmic content

- It forms the active-state basis `[U Hz Fz -HzFz]`, projects the relaxation-plus-pumping matrix into that four-state space, and displays the resulting matrix from Equation 2.

## Implementation structure

- A simulation of the matrix in Equation 2 of IK's paper on chemically
- amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011).
- Calculation time: seconds.
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)
- Relaxation theory
- Formalism and basis
- Spinach housekeeping
