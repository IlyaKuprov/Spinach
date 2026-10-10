# examples/kinetics/nonlinear/bimolecular_closures.m

Demonstrates `A+B -> C` with additive and product closures on two polarised one-proton reactants. Both closures use the reaction-record kernel and concentration-weighted initial states. RKMK4 endpoint states are compared with `ode45` at relative/absolute tolerances `1e-12`/`1e-14`; the plotted unweighted product coil detects the cross-reactant two-spin order, which is absent under additive closure.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.
