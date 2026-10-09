# tests/kernel/test_cwdm_flux.m

Checks intermolecular spin replacement as an explicit additive `A+B -> A+B` record. A two-spin molecule exchanges one spin with a one-spin pool; matching preserves the other molecular spin. Single-spin magnetisation transfers both ways, intramolecular orders involving the departing spin decay, and unaffected internal orders survive.

The test verifies zero concentration derivatives at independently sampled populations, magnetisation conservation, and agreement between the nonlinear production stepper and an exponential of the generator frozen at its invariant concentrations. This freezing is valid for time-independent additive replacement with no concentration-changing reactions; it is not a general mass-action approximation.

Detection and reference operator vectors explicitly use the `exact` method of the four-argument `coil_state` primitive.
