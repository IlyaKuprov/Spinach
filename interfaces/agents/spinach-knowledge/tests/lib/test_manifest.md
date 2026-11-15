# tests/lib/test_manifest.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

The regression registry maps stable kernel test IDs and descriptions to test functions. The existing CCE-family entries are retained, with kernel/clexpn_options and kernel/clexpn_transition added for the staged native input services. The new functions return standard result records and are executable through run_test and run_tests; they do not replace or weaken the existing CCE assertions.

The native retained-system registrations also include kernel/clexpn_subsystem and kernel/clexpn_descriptors. These exercise the compiler, observed-first physical mappings, cache invalidation, canonical operator/state descriptors and explicit populations. Existing registrations and assertions remain unchanged.

Native stage-2b registrations: kernel/clexpn_generators (generator signs, propagation and supported pairs) and kernel/clexpn_observation (actual reference, phase and detection).

## Native context and closed-system experiments

Registered kernel/clexpn_context and kernel/clexpn_experiments test the context, five-argument dispatch, closed-system experiments and recursive interlacing evaluator retention.
