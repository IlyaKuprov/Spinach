# examples/nmr_solids/cp_powder_static_nhh.m

- Signature: `cp_powder_static_nhh()`

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius sphere around it. Static powder simulation in a reduced (up to, and including four-spin correlations) Liouville space. Calculation time: minutes on a Tesla A100, much longer on CPU.

## Physical / mathematical content

This static 15N–1H cross-polarisation model has eight protons positioned around one 15N; the source describes the proton bath as scattered on a 2 Å-radius sphere. The simulation starts from the anisotropic equilibrium requested by `parameters.needs={'aniso_eq'}` at 298 K and observes the 15N response while both nuclei are irradiated in the doubly rotating frame.

## Numerical / algorithmic content

The calculation uses the sphten-liouv formalism with the IK-0 reduced basis, `inter_level=4`, and `sys.enable={'greedy'}` (the source comments that this example needs a GPU). The hard-contact pulse simulation is powder averaged on `rep_2ang_100pts_sph`; it uses 100 intervals of 10 μs and 50 kHz irradiation on both channels.

## Implementation structure

Builds the nine-spin system and fourth-level interaction-space basis, defines the two-channel RF and 15N detection operators, then calls `powder` with `cp_contact_hard` and plots the real 15N response versus cumulative contact time.
