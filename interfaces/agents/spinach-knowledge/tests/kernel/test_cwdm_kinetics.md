# tests/kernel/test_cwdm_kinetics.m

Tests direct-sum chemistry against analytic mass-action and spin-transport references (T8–T11). It constructs concentration-weighted initial states with `state` and `unit_state`, and detects product orders with `coil_state`.

The cases cover additive bimolecular arrival, arbitrary spin orders at independently sampled concentrations, zero concentration, reverse reactions, tracked spin-free sinks, repeated spin-free reactants, independent voxels, constant first-order exchange against a chemical rate-matrix exponential, time-dependent rates (including a profiled single evaluation shared by three unequal-population voxels), and product versus additive cross-reactant orders. RKMK4 convergence is measured against an `ode45` reference at relative/absolute tolerances `1e-12`/`1e-14`; halving the time step must reduce the error by a factor between 15 and 17, with finest-step error below `1e-10`.

A counted time-rate callback checks one evaluation for a three-voxel generator call; the resulting voxel derivatives must exactly match scaled single-voxel references.

This compact analytic network is not a substitute for regression of the migrated chemistry examples.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.

Zero numeric higher-order rates give a sparse zero generator and the same ordinary liquid/acquire FID as absent chemistry; a zero-valued rate callback still returns a dynamic handle.

The reporting regression combines column and row memberships in one reaction, checks the traced-spin text, and verifies exact equality with the ordinary-call generator.

Matched repeated products are rejected by identifier, while repeated spin-free products are checked against the corresponding mass-action stoichiometry.
