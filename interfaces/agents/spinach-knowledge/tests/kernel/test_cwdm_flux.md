# tests/kernel/test_cwdm_flux.m

Checks intermolecular spin replacement as an explicit additive `A+B -> A+B` record. A two-spin molecule exchanges one spin with a one-spin pool; matching preserves the other molecular spin. Single-spin magnetisation transfers both ways, intramolecular orders involving the departing spin decay, and unaffected internal orders survive.

The test verifies zero concentration derivatives at independently sampled populations, magnetisation conservation, and agreement between the nonlinear production stepper and an exponential of the generator frozen at its invariant concentrations. This freezing is valid for time-independent additive replacement with no concentration-changing reactions; it is not a general mass-action approximation.

Detection and reference operator vectors use the three-argument unweighted `coil_state` primitive.

Current interface note: flux rates/type live at inter.flux_rate and inter.flux_type, then spin_system.rlx.flux_rate and spin_system.rlx.flux_type. Selecting 'flux' in inter.relaxation activates the contribution; tests of prior automatic or chemistry-owned behavior are historical, not current API assertions.
