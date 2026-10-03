# tests/kernel/test_prop_cache.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_prop_cache.m)

`result=test_prop_cache()` checks shared-cache agreement with fresh propagation after changing cleanup, storage thresholds, and backend policy, including both policy insertion orders and worker retrieval. It uses the standard regression result helpers and an explicitly sized Spinach process pool. The backend partition check uses a small unscaled generator and does not require or validate GPU hardware.
