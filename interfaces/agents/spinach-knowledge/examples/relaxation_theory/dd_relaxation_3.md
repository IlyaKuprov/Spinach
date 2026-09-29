# examples/relaxation_theory/dd_relaxation_3.m

- MATLAB implementation: [examples/relaxation_theory/dd_relaxation_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/dd_relaxation_3.m)

- Signature: `dd_relaxation_3()`
- Run from MATLAB with no arguments; the function returns no MATLAB output and prints two rate ratios.
- Calculation time stated by the source: seconds.

## Purpose

This short extreme-narrowing comparison checks the dipolar longitudinal-relaxation rate ratio for a two-proton pair against a proton–deuteron pair. It tests the expected dependence on the square of the gyromagnetic-ratio ratio and on the partner-spin factor `S(S+1)`; it is a comparison, not a fitting or parameter-sweep routine.

## Model used

Both systems are evaluated at `sys.magnet = 14.1`, with the spins at `[0 0 0]` and `[0 0 2.50]`. The first system is `{'1H','1H'}`; the second is `{'1H','2H'}`. For each, the example sets Redfield relaxation, zero equilibrium, lab-frame retention, `inter.tau_c={1e-12}`, and a full `sphten-liouv` basis (`bas.approximation='none'`). The source gives no explicit unit annotation for the field, coordinates, or correlation time, so keep its values and Spinach unit conventions together rather than assigning units from this file alone.

## Calculation and interpretation

For each isotope pair the function builds the system and basis, constructs `R=relaxation(spin_system)`, forms the normalised longitudinal state `rho=state(spin_system,{'Lz'},{1})`, and obtains the scalar `rho'*R*rho`. It prints the ratio `R1pp/R1pd` and the textbook comparison

`[((1/2)(1/2+1))/(1(1+1))] * [spin('1H')/spin('2H')]^2`.

Agreement is the intended check of the extreme-narrowing scaling. The function only displays these values: it has no assertion, tolerance, returned rate, or automated pass/fail result. The comparison also relies on the particular two-spin geometry, relaxation settings, and basis above; changing them changes what is being tested.
