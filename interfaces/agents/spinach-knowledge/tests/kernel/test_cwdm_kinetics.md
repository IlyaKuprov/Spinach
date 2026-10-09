# tests/kernel/test_cwdm_kinetics.m

Tests direct-sum chemistry against analytic mass-action and spin-transport references (T8–T11). It constructs concentration-weighted initial states with `state` and `unit_state`, and detects product orders with `coil_state`.

The cases cover additive bimolecular arrival, arbitrary spin orders at independently sampled concentrations, zero concentration, reverse reactions, tracked spin-free sinks, repeated spin-free reactants, independent voxels, constant first-order exchange against a chemical rate-matrix exponential, time-dependent rates, and product versus additive cross-reactant orders. RKMK4 convergence is measured against an `ode45` reference at relative/absolute tolerances `1e-12`/`1e-14`; halving the time step must reduce the error by a factor between 15 and 17, with finest-step error below `1e-10`.

This compact analytic network is not a substitute for regression of the migrated chemistry examples.
