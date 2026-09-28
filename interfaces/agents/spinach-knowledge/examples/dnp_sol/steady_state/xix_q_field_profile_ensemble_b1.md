# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1.m

- Signature: `xix_q_field_profile_ensemble_b1()`

## Purpose

Simulates a steady-state XiX DNP field profile at Q band, averaged over an ensemble of electron Rabi frequencies. The source estimates a calculation time of minutes.

## Physical / mathematical content

The system contains an electron and a proton at a 3.5 Å separation in a 1.2142 T magnetic field, with a trityl electron g-tensor, an estimated proton chemical shift, and a temperature of 80 K. It uses `t1_t2` relaxation, including a proton longitudinal relaxation rate supplied by `r1n_dnp` that depends on orientation through `bet`. The detected observable is proton `Lz`.

## Numerical / algorithmic content

The script creates the spin system in an unrestricted `sphten-liouv` basis. It evaluates `xixdnp_steady` through `powder` on the `rep_2ang_800pts_sph` grid for 201 electron microwave offsets from −100 to 100 MHz. Each calculation uses 36 XiX blocks, 48 ns pulses, an inverted second-pulse phase, and the specified shot spacing. Five Gauss–Legendre points span electron nutation frequencies of 10–20 MHz; their weighted results are combined into an ensemble-averaged profile.

## Implementation structure

After configuring the spin system, relaxation, and experiment parameters, the function loops over the B1 quadrature points and runs the steady-state powder simulation at each point. It plots the real part of the averaged proton `Lz` expectation value against microwave resonance offset and saves `xix_q_field_profile_ensemble_b1.fig`.