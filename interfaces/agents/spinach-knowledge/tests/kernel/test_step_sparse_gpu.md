# tests/kernel/test_step_sparse_gpu.m

Verifies small-Hilbert `step` propagation using actual E8/E8 spin operators and independent angular-momentum rotation identities. Real and complex Hermitian generators are tested with sparse/full CPU inputs, GPU promotion, existing GPU inputs, numeric and cell states, and zero duration. CPU comparisons remain active without a GPU; unavailable GPU checks are explicitly reported as skipped. The function returns a standard regression result and is registered as `kernel/step_sparse_gpu`.
