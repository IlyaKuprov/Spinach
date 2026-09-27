# examples/spin_chemistry/cidnp_pumping_2.m

- Signature: `cidnp_pumping_2()`

## Purpose

A simulation of Figure 2A in IK's paper on chemically amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011). Calculation time: seconds.

## Physical / mathematical content


## Physical / mathematical content

- At 14.1 T, the model is a 1H-19F pair with zero isotropic shifts, 50 Hz scalar coupling, and a fluorine CSA tensor with principal values `[-47 -16 63]`; the spins are separated by `2.60` along y.
- Redfield relaxation retains secular terms, uses `tau_c=110e-12`, and is configured for equilibrium magnetisation at `298 K`. Pumping terms act on the proton and fluorine z-magnetisations; the proton also receives an additional relaxation term with coefficient `3.0`.

## Numerical / algorithmic content

- Starting from `unit_state + 2*Hz + 2*Fz`, the script evolves for 40 steps of `0.1 s` using three detection channels, `Fz`, `-HzFz`, and `Hz`. It plots these signals over 4 s.

## Implementation structure

- A simulation of Figure 2A in IK's paper on chemically amplified
- NOEs (https://doi.org/10.1016/j.jmr.2004.01.011).
- Calculation time: seconds.
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)
- Relaxation theories
- Formalism and basis
- Spinach housekeeping
