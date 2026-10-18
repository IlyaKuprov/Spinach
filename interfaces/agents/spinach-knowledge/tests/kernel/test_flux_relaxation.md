# tests/kernel/test_flux_relaxation.m

`result=test_flux_relaxation()` checks the explicit flux relaxation contract with a two-spin system. It verifies an independently specified directed-transfer column, enabled and disabled operation, zero rates, top-level input and internal `rlx` ownership, summary gating, and separation from kinetics. Nonzero stored flux is also tested while disabled in Zeeman Liouville space.

The combined-mechanism cases check additivity with T1/T2 relaxation and with a nonzero anisotropic Redfield generator. Comparisons use exact zeros or tolerances scaled to the generator norm and floating-point precision. Incomplete flux parameter pairs and retired `inter.chem.flux_*` inputs must be rejected. Merge checks preserve the top-level fields and reject incorrect constituent dimensions; spin deletion must remove the corresponding internal rate row and column.

`test_cwdm_flux` separately covers both correlation-transfer models against an independently assembled full generator, retention and thermalisation, SRSK recursion, cross-substance reduction, and zero-rate behaviour. `test_cwdm_contexts` checks that imaging only uses caller-supplied relaxation operators.
