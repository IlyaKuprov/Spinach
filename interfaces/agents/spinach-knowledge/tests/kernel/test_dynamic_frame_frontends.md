# tests/kernel/test_dynamic_frame_frontends.m

- Signature: `result=test_dynamic_frame_frontends()`

## Purpose

Tests algebraic limits of the `carrier()`, `frqoffset()`, `rotframe()`, `average()`, and `orientation()` front ends on compact systems.

## Physical / mathematical content

- In a heteronuclear `1H`/`13C` spherical-tensor Liouville-space system, checks that the `carrier()` commutator and anticommutator forms equal the difference and sum of its left and right forms; that the all-spin commutator equals the sum of isotope-specific carriers; and that the proton carrier equals its free-particle Larmor frequency times the `Lz` commutator operator.
- Checks `frqoffset()` against `2*pi*25.0*Lz(1H) - 2*pi*40.0*Lz(13C)`, both for independent channels and for a duplicated `1H` channel with the same offset.
- Checks that zeroth-order `rotframe()` removes the carrier from a one-spin Hilbert-space Hamiltonian, leaving the transverse perturbation `H1`.
- Checks that `average()` leaves an unmodulated Hamiltonian `H0 = [0 1; 1 0]` unchanged for first-, second-, and third-order average-Hamiltonian (`ah_*`) and Krylov-Bogolyubov (`kb_*`) branches, and for `matrix_log`, with `omega = 2*pi*1000`.
- Checks `orientation()` at non-zero Euler angles `[0.2 0.3 -0.4]` against an explicit rank-one and rank-two rotational-basis contraction: `D1(1,3)*Q{1}{1,3} + D1(3,1)*Q{1}{3,1} + D2(3,3)*Q{2}{3,3}`, Hermitian-symmetrized, where `D1` and `D2` come from `wigner()` at those angles.

## Numerical / algorithmic content

- The `carrier()` checks use relative tolerance `1e-8` and absolute tolerance `1e-14`; `frqoffset()` uses `1e-12` for both; `rotframe()` uses `1e-10` and `1e-14`; `average()` uses `1e-10` and `1e-12`; and `orientation()` uses `1e-14` for both.
- The common `1H`/`13C` system uses a `14.1` magnet, scalar Zeeman values `1.0` and `2.0`, a scalar coupling of `10.0`, and an untruncated `sphten-liouv` basis. The rotating-frame check instead uses a one-proton `zeeman-hilb` system with zero scalar Zeeman shift and a Hermitian transverse perturbation formed from `1e3*Lx`.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Creates a test result for `kernel/dynamic_frame_frontends`, then runs separate algebraic checks for `carrier()`, `frqoffset()`, `rotframe()`, `average()`, and `orientation()`.

ilya.kuprov@weizmann.ac.il