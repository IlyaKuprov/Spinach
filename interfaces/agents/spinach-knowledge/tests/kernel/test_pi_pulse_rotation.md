# tests/kernel/test_pi_pulse_rotation.m

- Signature: `result=test_pi_pulse_rotation()`

## Purpose

Tests that a hard `pi` rotation about `X` inverts longitudinal magnetisation: `Lz -> -Lz`.

## Physical / mathematical content

An active rotation by `pi` around `X` maps `z` to `-z`.

## Numerical / algorithmic content

- Builds a one-proton Hilbert-space spin system using the `zeeman-hilb` formalism with no approximation, zero magnet setting, and zero scalar Zeeman interaction.
- Constructs `Lx` and `Lz`, then computes `rho_obs=step(spin_system,Lx,Lz,pi)`.
- Compares `rho_obs` with `-Lz` using absolute and relative tolerances of `1e-14`.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announces the hard pi-pulse rotation test and creates its test result.
- Builds the spin system, applies the X rotation, and checks that Lz is inverted.