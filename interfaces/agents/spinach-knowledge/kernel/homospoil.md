# kernel/homospoil.m

- Signature: `rho=homospoil(spin_system,rho,zqc_flag)`

## Purpose

Emulates a strong homospoil pulse: only states with zero frequency relative to the carrier frequencies survive; chemical shifts are not considered. In `sphten-liouv`, `zqc_flag` controls whether zero-quantum coherences are retained (`keep`) or destroyed (`destroy`). In `zeeman-liouv` and `zeeman-hilb`, zero-quantum coherences are always destroyed, leaving only the density-matrix diagonal.

## Parameters / inputs

- `rho`: a state vector or a horizontal stack of state vectors.
- `zqc_flag`: `'keep'` or `'destroy'`.

## Output

- `rho`: state vector(s) containing the retained longitudinal states and, when requested, zero-quantum coherences.

## Formalism support

Fokker-Planck direct products are supported in the Liouville-space formalisms.

## Source documentation

https://spindynamics.org/wiki/index.php?title=homospoil.m
