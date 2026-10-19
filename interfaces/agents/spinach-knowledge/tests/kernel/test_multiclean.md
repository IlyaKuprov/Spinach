# test_multiclean.m

Checks memory-bounded matrix multiplication against native multiplication followed by representation-preserving `clean_up`. It covers sparse/full combinations, real and complex values, single and double precision, scalar and empty shapes, zero-tolerance bypass, and rounding boundaries. Dyadic examples isolate cleanup semantics from legitimate blocked-BLAS summation differences.

`test_multiclean_gpu` applies the corresponding checks to supported gpuArray products; unavailable GPU coverage is reported as SKIP. `test_cleanup_storage` verifies that the helper’s global storage conversion does not round a fused result twice and preserves custom density thresholds and bypass semantics.

The CPU and GPU suites cover NaN and disabled-cleanup native arithmetic, custom density and size thresholds, and whole-result representation rather than panel-local density. CPU tests also exercise the unfused polyadic decoupling fallback.
