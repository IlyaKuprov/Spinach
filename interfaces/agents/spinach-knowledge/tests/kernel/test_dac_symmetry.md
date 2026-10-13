# tests/kernel/test_dac_symmetry.m

Tests the shared Diels–Alder molecular input with its shipped IK-2 bases and unchanged interactions. The solvent-only S3 declaration uses substance-local indices. Each proposed S2 group (global spins 5–6 and 25–26) is independently passed through `basis` and must be rejected because its full Zeeman tensors are unequal.

For the accepted solvent symmetry, the test compares basis descriptors to a no-symmetry build, verifies orthonormality and closure of the fully symmetric projector, and compares complex FIDs using projected preparation, detection, Hamiltonian, and the original solvent relaxation. Reacting spins are traced out only for the inexpensive solvent relaxation calculation. The FID comparison covers the solvent block, not the complete reacting-mixture spectrum or the hydrodynamic examples.
