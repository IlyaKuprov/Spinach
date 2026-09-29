# examples/singlet_states/dipolar_singlet.m

- Signature: `dipolar_singlet()`

## Purpose

Tests the action of a liquid-state dipolar Redfield relaxation superoperator on a two-proton singlet. The source describes the singlet as immune to dipolar relaxation, but the program itself prints a norm rather than asserting a tolerance or reporting a measured lifetime.

## Physical model

The two `1H` spins are placed at coordinates `[0.0 0.0 0.0]` and `[0.5 0.6 0.7]` in a `14.1 T` field. The source does not label coordinate units. It requests Redfield relaxation with zero equilibrium, lab-frame relaxation terms, and a correlation time of `5e-9 s`. Integration and zero tolerances are each `1e-5`; the proximity cutoff is `4.0` (no unit is specified).

## Calculation and observable

The calculation uses the unapproximated `sphten-liouv` basis, constructs the relaxation superoperator `R`, normalises the singlet state of spins 1 and 2, and displays `norm(R*S)`. This is a relaxation-action norm for that normalised state, not a simulated time trace or a measured singlet lifetime. The source gives no numeric output value.

## Source

[examples/singlet_states/dipolar_singlet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/dipolar_singlet.m)
