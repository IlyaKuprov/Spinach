# examples/relaxation_theory/csa_antisymm_1.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/csa_antisymm_1.m)

This example estimates longitudinal and transverse relaxation rates for one 13C nucleus from a Redfield relaxation matrix and compares those projections with the textbook CSA routine rlx_csa. It is a numerical model comparison, not an experimental benchmark. The field is 14.1 T; the chemical-shielding matrix (ppm) is [100 20 15; 20 0 30; 25 10 -30], which is nonsymmetric and therefore includes an antisymmetric component. The correlation time is 50 ps (50e-12 s).

The source requests Redfield relaxation with zero equilibrium and lab-frame retention, using sphten-liouv without approximation. It computes R=relaxation(spin_system), forms longitudinal and transverse estimates by projecting R onto Lz and L+ respectively, and applies a minus sign to each normalised projection. The same field, isotope, shielding matrix and correlation time are passed to rlx_csa for the textbook comparison. It prints both sets of rates; the source does not supply a measured value or state a numerical agreement claim. No cross-correlation term is explicitly selected.
