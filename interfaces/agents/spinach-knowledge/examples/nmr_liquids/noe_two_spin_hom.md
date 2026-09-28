# examples/nmr_liquids/noe_two_spin_hom.m

- Signature: `noe_two_spin_hom()`

## Purpose

Nuclear overhauser effect in a homonuclear two-spin system in the long correlation time case. Calculation time: seconds.

## Physical / mathematical content

- The model is a pair of proton spins separated by 2.00 Å, with zero isotropic Zeeman offsets. Redfield relaxation uses a 1 ns correlation time, 298 K, and the Di Bari equilibrium convention.
- One proton's longitudinal magnetization is inverted from thermal equilibrium. The two detected longitudinal components show the relaxation-mediated NOE response in the long-correlation-time regime.

## Numerical / algorithmic content

- The full spherical-tensor Liouville basis is used without a basis approximation; the Redfield superoperator retains the kite-selected terms.
- Multichannel relaxation evolution is sampled at 0.01 s intervals for 1000 intervals, spanning 0–10 s.

## Implementation structure

- Build the two-proton system at 14.1 T, form its basis and Redfield relaxation superoperator, then calculate equilibrium.
- Invert spin 1, propagate with both proton `Lz` operators as observation channels, and plot the longitudinal magnetizations labelled Proton A and Proton B.
