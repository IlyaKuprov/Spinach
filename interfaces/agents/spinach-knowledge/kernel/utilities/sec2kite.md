# kernel/utilities/sec2kite.m

- Signature: `R=sec2kite(spin_system,R)`

## Purpose

Converts a secular relaxation superoperator into Redfield kite form by dropping all non-longitudinal cross-relaxation processes. Useful when the relaxation superoperator is huge but TROSY-like effects are negligible.

## Parameters / inputs

- `spin_system` — spin system; its basis is used to identify longitudinal product states.
- `R` — relaxation superoperator.

## Outputs

- `R` — relaxation superoperator retaining self-relaxation and longitudinal cross-relaxation terms.

## Method

The function identifies longitudinal product states from `spin_system.bas.basis`, then retains matrix entries whose row and column both correspond to longitudinal product states, as well as all diagonal entries. Other entries are set to zero. It reports the number of nonzero entries before and after conversion.

## Requirements and caveats

- Requires the `sphten-liouv` formalism and a numeric, square `R`.
- Cannot proceed if `norm(R*unit_state(spin_system),2)>1e-10`; this is treated as an indication that `R` has been thermalised.
- Non-longitudinal cross-relaxation processes are discarded, so use this conversion when TROSY-like effects are negligible.

## Contact and reference

- ilya.kuprov@weizmann.ac.il
- <https://spindynamics.org/wiki/index.php?title=sec2kite.m>