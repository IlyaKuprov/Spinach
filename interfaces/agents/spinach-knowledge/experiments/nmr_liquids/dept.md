# experiments/nmr_liquids/dept.m

- Signature: `fid=dept(spin_system,parameters,H,R,K)`

## Purpose

DEPT pulse sequence; see [DOI 10.1016/0022-2364(82)90286-4](https://doi.org/10.1016/0022-2364(82)90286-4).

## Physical / mathematical content

- The sequence begins from isotropic thermal equilibrium, applies a 90-degree pulse to spin 2, and uses three J-coupling evolution periods of `abs(1/(2*parameters.J))`. The source constructs two phase alternatives with x/y pulses, combines them, then applies the `parameters.beta` editing pulse before the final coupling period.
- It decouples spin 2 and detects the FID on spin 1. The effective Liouvillian is `L = H + 1i*R + 1i*K`.

## Numerical / algorithmic content

- A single acquisition dimension uses dwell time `1/parameters.sweep` and `parameters.npoints` points. The implementation requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=dept(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: sweep width in Hz.
- `parameters.npoints`: number of acquisition points.
- `parameters.spins`: two spin labels `{F1,F2}`, e.g. `'13C'` and `'1H'`.
- `parameters.J`: working J-coupling in Hz.
- `parameters.beta`: selection-pulse angle in radians.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices received from the context function.

## Outputs

- `fid`: free induction decay.
- DEPT135 gives CH and CH3 signals opposite in phase to CH2; DEPT90 gives only CH signals; DEPT45 gives positive CH, CH2, and CH3 signals. Quaternary carbons do not appear.
- Use isotope dilution to generate carbon isotopomers; see `dilute.m`.

[Spinach Wiki: dept.m](https://spindynamics.org/wiki/index.php?title=dept.m)
