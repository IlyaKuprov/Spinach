# examples/fundamentals/state_tests/thermal_equilibrium_3.m

- Signature: `thermal_equilibrium_3()`
- Source: [`examples/fundamentals/state_tests/thermal_equilibrium_3.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/thermal_equilibrium_3.m)

## Purpose

Checks Spinach's isotropic thermal-equilibrium state against textbook Boltzmann-population expressions for the three spins, across three basis formalisms.

## Model and setup

The example uses a three-spin trityl/proton system: one electron and two `1H` spins. It sets `sys.magnet=0.34` (described in the source as an X-band magnet), anisotropic Zeeman principal values `[2.00319 2.00319 2.00258]` for the electron and proton shift guesses `[0 0 5]` and `[0 5 0]`. The Euler-angle vectors are `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees (the source scales them by `pi/180`); the three Cartesian coordinate vectors are `[0 0 0]`, `[0 3.5 0]`, and `[2.475 2.475 0]`. The source sets the spin-temperature parameter to `80` without spelling out its unit.

## Use and comparison

The basis uses `approximation='none'` and successively tests `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`. For each, it forms `equilibrium(spin_system)`, evaluates each spin's longitudinal expectation value (`Lz`),, and compares those values with expressions built from `levelpop` for `E` and `1H`. The stated relative-error limit is `1e-3`; a larger discrepancy triggers an error.

## Output and limits

This is an assertion-style test, not a spectrum simulation. The source contains a success message, but that message is not evidence that the test was run here. The field, temperature, and coordinate values above are reported as source inputs; the file does not declare units for the temperature or coordinates.
