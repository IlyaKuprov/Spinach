# examples/spin_chemistry/cidnp_flash_acquire.m

- Signature: `cidnp_flash_acquire()`

## Purpose

A model of the CIDNP magnetisation pumping process described in IK's paper: The system is pumped for 0.5 seconds, and then allowed to relax. Calculation time: seconds. Miguel Mompean Ilya Kuprov

## Physical / mathematical content

- At 14.1 T, the model contains 1H and 19F with zero isotropic shifts, a 50 Hz scalar coupling, and a fluorine CSA tensor with principal values `[-47 -16 63]`. The spins are separated by `2.60` along the y coordinate.
- Relaxation is Redfield with secular terms retained, zero equilibrium, and correlation time `110e-12`. The relaxation matrix is thermalised to the sum of the proton and fluorine z-magnetisations; light-induced pumping terms are added for each spin.

## Numerical / algorithmic content

- Starting from `unit_state + Hz + Fz`, the script propagates under illumination for 50 steps of `0.01 s`, then without the light-pumping terms for 500 further steps of `0.01 s`. It plots the time courses of `Fz`, `Hz`, and `HzFz` over 5.5 s.

## Implementation structure

- A model of the CIDNP magnetisation pumping process described in
- IK's paper:
- The system is pumped for 0.5 seconds, and then allowed to relax.
- Calculation time: seconds.
- Miguel Mompean
- Ilya Kuprov
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)
