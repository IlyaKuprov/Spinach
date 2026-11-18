# tests/kernel/test_clexpn_orientation.m

Registered test kernel/clexpn_orientation covers deferred context and powder dispatch without a full-bath basis, normalized explicit quadratures, both noncommuting averaging orders, raw-reference normalization, combined sampled-state/orientation reduction and standard errors, integer-weight normalization, shipped named-grid lookup, zero-weight omission, malformed descriptor and grid rejection, assumptions agreement, and the unchanged ordinary five-argument pulse-sequence path.

Additional checks compare every orientation-resolved amplitude, error bound, and validity mask against direct clexpn calls; compare serial and cluster-parallel composition bitwise; run composition inside an existing worker with serial inner scheduling; and verify that rotation changes an anisotropic Zeeman pump selection and the corresponding independently evaluated responses.
