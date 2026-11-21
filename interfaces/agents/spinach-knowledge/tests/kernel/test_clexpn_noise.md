# test_clexpn_noise

Syntax: `result=test_clexpn_noise()`

Registered native contract tests for unequal conditional means, ensemble connected subtraction, quantum ordering, raw first moments, and empty-reference inclusion-exclusion.

The noise example uses native operator descriptors. The cosine reference and original transverse operator are retained in the example-test capture. Native noise tests compare descriptors and explicit Hermitian matrices before connected ensemble reduction and reject invalid descriptor fields.

Direct registered-evaluator tests reject noise input paired with a foreign experiment schedule. The exact operator 1e-20*[0 1;0 0] is rejected using the operator-relative Hermiticity contract, independently of physical units.

The sampled unequal-mean fixture also checks the centred-contribution standard error exactly against the empirical covariance expression; it is not the uncertainty of uncentred second moments.

A polarised, relaxing spin with a complex transverse fluctuation checks both quantum orderings against exponential population and coherence decay. Connected and raw cases distinguish the first-moment cost, and rescaling the fluctuation tests squared-unit error scaling.
