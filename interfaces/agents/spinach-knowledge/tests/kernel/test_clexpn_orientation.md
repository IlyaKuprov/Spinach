# tests/kernel/test_clexpn_orientation.m

Registered test `kernel/clexpn_orientation` covers native single-orientation and grid calls, weighted and per-orientation outputs, both averaging orders, shared sampled spin states, preservation of bounds and masks, and checkpoint restrictions.

Additional checks compare every orientation-resolved amplitude, error bound, and validity mask against direct clexpn calls; compare serial and cluster-parallel composition bitwise; run composition inside an existing worker with serial inner scheduling; and verify that rotation changes an anisotropic Zeeman pump selection and the corresponding independently evaluated responses.
