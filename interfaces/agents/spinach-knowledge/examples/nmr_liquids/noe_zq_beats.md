# examples/nmr_liquids/noe_zq_beats.m

- Signature: `noe_zq_beats()`

## Purpose

A relaxation-time simulation of zero-quantum beats in the Overhauser effect for a strongly coupled homonuclear two-spin system. It follows longitudinal signals after one spin is inverted; it does not run a NOESY pulse sequence. The source estimates a calculation time of seconds.

## Spin system and Liouvillian

The two spins are `1H`. Their isotropic Zeeman scalars are 0.0 and 0.01 ppm; the scalar-coupling matrix contains 3.0 Hz for the pair. Coordinates are `[0 0 0]` and `[0 0 2]` Angstrom. The field is set by `sys.magnet=14.1`. The basis is spherical-tensor Liouville with no approximation. Redfield relaxation uses the `dibari` equilibrium convention, `secular` retention, temperature 298, and `tau_c={1e-9}`; the proximity cut-off is 4.0. The propagated Liouvillian is the assumed NMR Hamiltonian plus `1i*relaxation(spin_system)`.

## Preparation and observables

Thermal equilibrium is calculated using the lab-frame Hamiltonian. The initial state inverts spin 1's `Lz` component relative to equilibrium. Multichannel propagation detects both spins' `Lz` operators, with step parameter `1e-2`, 1000 steps, and `multichannel` mode. The code plots the real part of the answer against `linspace(0,10,1001)`, with time labelled in seconds.

## Source

[examples/nmr_liquids/noe_zq_beats.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noe_zq_beats.m)
