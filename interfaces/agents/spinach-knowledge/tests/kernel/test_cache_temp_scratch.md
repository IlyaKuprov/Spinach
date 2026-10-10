# tests/kernel/test_cache_temp_scratch.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cache_temp_scratch.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cache_temp_scratch.m)

## Purpose

Regression test for cache management in a temporary scratch directory. It checks `cacheman()` and `wipe_cache()` against a temporary scratch directory, loads shipped small cache tables, exercises the on-demand SLE operator cache, and smoke-tests the read-only `sniff()` integrity pass. The stated requirement is that `cacheman()` and `wipe_cache()` must only affect Spinach cache records in the configured scratch directory.

## What is tested

The cache test isolates a Spinach scratch cache from an unrelated ordinary file in a temporary directory. At the 365-day retention threshold, `cacheman()` must keep a fresh Spinach cache entry and leave the ordinary file untouched. `wipe_cache()` removes the Spinach entry but not the ordinary file; with a zero-day threshold, `cacheman()` removes an expired Spinach cache directory. These checks distinguish retention and selective deletion from indiscriminate scratch cleanup.

The shipped two-level `st_product_table`, `ist_product_table`, and `bos_product_table` records must remain loadable, with left and right arrays of size `[4 4 4]`. The cached single-rank SLE operators have `10×10` Cartesian matrices, a `[10 3]` Wigner basis and a `[1 2]` `D` cell (`D{1}` empty, `D{2}` of size `[5 5]`). A read-only `sniff('none')` smoke check must emit an integrity report.

## Inputs and outputs

- Syntax: `result=test_cache_temp_scratch()`
- Output: `result` — regression test result with explanatory messages.
- No inputs.

## References

- `cacheman()` — cache retention management exercised against the scratch directory.
- `wipe_cache()` — cache wipe routine checked for selective deletion.
- `st_product_table()`, `ist_product_table()`, `bos_product_table()` — shipped product table caches loaded with level argument `2`.
- `sle_operators()` — on-demand SLE operator cache exercised with arguments `(1,2)`.
- `sniff()` — read-only integrity pass smoke-tested with `'none'`.
- `new_test_result()`, `test_true()` — test harness helpers.
