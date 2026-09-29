# examples/nmr_liquids/noe_two_spin_het.m

- Signature: `noe_two_spin_het()`

## Purpose

A relaxation-time trajectory illustrating the short-correlation-time nuclear Overhauser effect in a heteronuclear two-spin system. This is not a NOESY pulse-sequence simulation: the example prepares an inverted longitudinal state and follows both spins' longitudinal magnetisation under Redfield relaxation. The source estimates a calculation time of seconds.

## Spin system and relaxation model

The isotopes are `1H` and `13C`; both isotropic Zeeman scalars are set to zero. Their coordinates are `[0 0 0]` and `[0 0 1.03]`; the source gives no coordinate unit. The field is set by `sys.magnet=14.1`. The basis uses spherical-tensor Liouville formalism with no approximation. Redfield relaxation uses the `dibari` equilibrium convention, `kite` retention, temperature 298, and `tau_c={100e-12}`.

## Preparation and observable channels

After constructing the relaxation superoperator and thermal equilibrium state, the code inverts spin 1's `Lz` component relative to equilibrium. Multichannel evolution uses the `Lz` operators for spins 1 and 2 as separate observation channels. The call uses `1i*R`, step parameter `1e-2`, 400 steps, and `multichannel` mode; the plotted time coordinate is `linspace(0,4,401)` and is labelled in seconds. The plot labels the traces Proton and Carbon.

## Source

[examples/nmr_liquids/noe_two_spin_het.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noe_two_spin_het.m)
