# test_multiclean.m

Compares panel multiplication with native multiplication followed by `clean_up`, allowing numerical roundoff rather than requiring bitwise or storage identity. Coverage includes sparse/full combinations, real and complex values, single and double precision, non-dyadic products, scalar and empty shapes, custom operators, and cleanup policy switches. Nonfinite and non-floating tolerances must be rejected.

`test_multiclean_gpu` checks supported GPU products and mixed CPU/GPU inputs; unavailable GPU coverage is reported as skipped. `test_cleanup_storage` checks ordinary cleanup storage policy and numerical agreement. These numerical tests do not establish a hard memory bound or a workload-specific peak-memory reduction.
