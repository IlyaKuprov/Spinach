# tests/kernel/test_spin_table.m

Source: [tests/kernel/test_spin_table.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_spin_table.m)

`test_spin_table()` protects the physical-data boundary of `spin`: an unknown property must not become a guessed simulation input. It rebuilds an isolated copy of the literature TSV and compares the complete table, including qualifiers and metadata, with the shipped MAT payload.

The test distinguishes spin-zero nuclei, stable tantalum's metastable state, abundance intervals, unsigned moments, unknown labels, and tentative nuclear assignments. Particle checks cover charge/CPT signs, mean-life to half-life conversion, and the spin-three-halves Omega moment conversion, and its unknown rather than fabricated-zero quadrupole moment. The isotope cutoff excludes subsecond values, estimates, and upper limits, retains unresolved lower limits and the exact one-second boundary, and rejects removed keys including neutral-atom `229Th_m`; longer alternative crystal lifetimes do not replace its adopted neutral-atom lifetime. Abstract modes have typed empty metadata, and editing returned metadata must not mutate the persistent table.

A direct Hamiltonian comparison verifies that a spin-zero spectator does not alter the proton's Hilbert-space dynamics. The test is registered as `kernel/spin_table`; it does not certify the completeness or bibliographic accuracy of the literature inventory.
