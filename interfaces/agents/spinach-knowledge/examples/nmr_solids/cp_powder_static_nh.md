# examples/nmr_solids/cp_powder_static_nh.m

- Signature: `cp_powder_static_nh()`

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame. Static powder simulation. Calculation time: seconds

## Physical / mathematical content

The two-spin system contains 15N and 1H at the specified 1.05 Å separation; both isotropic Zeeman shifts are set to zero. The experiment applies simultaneous RF irradiation to the two nuclei for cross-polarisation, and the static powder average is detected on 15N.

## Numerical / algorithmic content

The example uses the full sphten-liouv basis (approximation='none') and powder with cp_contact_hard. It averages over rep_2ang_6400pts_sph, with 100 contact intervals of 10 μs each; both irradiation channels are set to 5×10⁴. The plotted signal is the real part of the 15N FID against cumulative contact time.

## Implementation structure

Defines the spin system and basis, supplies the RF operators, powers, detection state, powder grid, time steps, and aniso_eq requirement, then runs the static powder simulation and plots the 15N response.
