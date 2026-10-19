# test_multiclean.m

Checks memory-bounded matrix multiplication against native multiplication followed by representation-preserving `clean_up`. It covers sparse/full combinations, real and complex values, single and double precision, scalar and empty shapes, zero-tolerance bypass, and rounding boundaries. Dyadic examples isolate cleanup semantics from legitimate blocked-BLAS summation differences.

`test_multiclean_gpu` applies the corresponding checks to supported gpuArray products; a usable GPU is required. `test_cleanup_storage` verifies that caller-level storage conversion does not round a fused result twice and preserves custom density thresholds and bypass semantics.
