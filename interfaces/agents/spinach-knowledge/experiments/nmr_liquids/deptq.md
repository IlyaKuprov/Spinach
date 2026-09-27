# experiments/nmr_liquids/deptq.m

- Signature: `fid=deptq(spin_system,parameters,H,R,K)`

## Purpose

DEPTQ pulse sequence; see [DOI 10.1006/jmre.1998.1595](https://doi.org/10.1006/jmre.1998.1595). This is the DEPTQ135-style variant with a fixed first proton pulse; `parameters.beta` controls the final proton editing pulse.

## Physical / mathematical content

- The source starts from isotropic thermal equilibrium, applies a carbon pulse, and uses J-coupling intervals of `abs(1/(2*parameters.J))` with phase-alternated carbon and proton pulses. The final beta pulse and detection pulse precede proton decoupling and FID acquisition on the carbon spin.
- Unlike `dept.m`, this sequence allows quaternary-carbon signals. The effective Liouvillian is `L = H + 1i*R + 1i*K`.

## Numerical / algorithmic content

- The single-dimension dwell time is `1/parameters.sweep`; acquisition uses `parameters.npoints` points. The implementation requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=deptq(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: sweep width in Hz.
- `parameters.npoints`: number of acquisition points.
- `parameters.spins`: two spin labels `{F1,F2}`, e.g. `'13C'` and `'1H'`.
- `parameters.J`: working J-coupling in Hz.
- `parameters.beta`: angle of the selection pulse, in radians.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices received from the context function.

## Outputs

- `fid`: free induction decay.
- Use isotope dilution to generate carbon isotopomers; see `dilute.m`.

[Spinach Wiki: deptq.m](https://spindynamics.org/wiki/index.php?title=deptq.m)
