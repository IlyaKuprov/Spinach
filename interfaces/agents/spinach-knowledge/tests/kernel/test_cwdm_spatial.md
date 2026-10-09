# tests/kernel/test_cwdm_spatial.m

The two-cell test uses two additive bimolecular channels and a spin-free solvent. Removing all spins must retain five scalar chemical blocks. Independent mass-action derivatives are compared with both concentration-only and full-spin kernel generators, including transport between cells and a nonuniform solvent population. Full-spin propagation must give the same unit-coordinate trajectory as the traced model to a whole-vector absolute tolerance of `1e-12`.

The test also distinguishes the original reacting-flow frozen concentration generator from the new additive generator. Both give the same derivative at the assembly state, but their product-unit source allocations differ: the original assigns the source to one reactant, while the kernel shares it equally. Consequently their finite frozen exponentials are not identical. This is a compact integration test, not full-chip T20 acceptance or a WP0 comparison.

The homogeneous two-channel case also compares frozen stock and kernel allocations with the exact bimolecular concentration solution on successively halved steps. Both allocations converge to the same mass-action trajectory; the test separately asserts the plan-prescribed equal unit-source sharing of the additive kernel. Finite frozen steps are not claimed to reproduce the exact ODE.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.
