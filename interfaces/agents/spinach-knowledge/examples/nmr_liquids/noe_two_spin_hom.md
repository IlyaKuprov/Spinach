# examples/nmr_liquids/noe_two_spin_hom.m

- Signature: `noe_two_spin_hom()`

## Purpose

A relaxation-time trajectory for the nuclear Overhauser effect in a homonuclear two-spin system in the long-correlation-time case. One longitudinal spin component is inverted and both spins' longitudinal magnetisations are followed; this is not a NOESY acquisition. The source estimates a calculation time of seconds.

## Spin system and relaxation model

Both isotopes are `1H`, both isotropic Zeeman scalars are zero, and the source assigns coordinates `[0 0 0]` and `[0 0 2.00]` without giving a coordinate unit. It sets `sys.magnet=14.1`. The basis uses spherical-tensor Liouville formalism with no approximation. Redfield relaxation uses the `dibari` equilibrium convention, `kite` retention, temperature 298, and `tau_c={1e-9}`.

## Preparation and observable channels

The code constructs the relaxation superoperator and thermal equilibrium state, then inverts spin 1's `Lz` component relative to equilibrium. Multichannel evolution under `1i*R` detects the two spins separately with their respective `Lz` operators, using step parameter `1e-2` and 1000 steps. The trajectory is plotted against `linspace(0,10,1001)`, labelled in seconds, with the traces identified as Proton A and Proton B.

## Source

[examples/nmr_liquids/noe_two_spin_hom.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noe_two_spin_hom.m)
