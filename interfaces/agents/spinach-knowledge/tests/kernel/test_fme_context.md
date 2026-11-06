# tests/kernel/test_fme_context.m

Checks the new high-frequency calculation against existing Spinach Fokker-Planck machinery, without replacing the context or its generator extractor. A two-proton secular dipolar MAS system is constructed in full spherical-tensor Liouville space. `drifts(spin_system,@singlerot,...)` supplies the native enlarged generator.

The five spin-space MAS harmonics come from the separate native `rotor_stack` helper. For an odd rotor grid, the spectral derivative has zero diagonal, so reconstruction against the Fokker-Planck diagonal phase blocks checks the two native generation paths before comparing dynamics. A localised initial zero-phase state and a complex transverse detector make the propagation comparison phase sensitive. Both endpoint kicks are retained, using the signed native frequency `-2*pi*parameters.rate`.

Native references at rotor harmonic cutoffs eight and ten must agree independently of the effective-generator comparison. This is a specific secular rank-two MAS integration test, not a general extractor for arbitrary flow, diffusion, or nonuniform spatial grids.
