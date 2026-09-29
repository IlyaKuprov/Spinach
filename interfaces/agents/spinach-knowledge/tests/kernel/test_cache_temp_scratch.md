# tests/kernel/test_cache_temp_scratch.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cache_temp_scratch.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cache_temp_scratch.m)

## Purpose

Regression test for cache management in a temporary scratch directory. It checks `cacheman()` and `wipe_cache()` against a temporary scratch directory, loads shipped small cache tables, exercises the on-demand SLE operator cache, and smoke-tests the read-only `sniff()` integrity pass. The stated requirement is that `cacheman()` and `wipe_cache()` must only affect Spinach cache records in the configured scratch directory.

## Behaviour

- Announces the test target with `fprintf('TESTING: Temporary scratch cache management\n')` and initialises a test result via `new_test_result('kernel/cache_temp_scratch', ...)`.
- Creates an isolated scratch directory with `tempname(tempdir)` and `mkdir`, registering `onCleanup(@()local_remove_dir(scratch))` for cleanup.
- Builds a minimal spin system with `sys.output='hush'`, `sys.scratch=scratch`, empty `sys.enable` and `sys.disable` cell arrays, and `tols.cache_mem=365`.
- Ensures `cacheman()` uses a small process pool instead of auto-starting a large one: if `gcp('nocreate')` returns empty, it calls `parpool('Processes',1)`.
- Saves two files into the scratch directory: `spinach_fresh.mat` and `ordinary_file.mat`, each containing `payload=1`.
- With the long retention horizon (`cache_mem=365`), verifies `cacheman()` keeps the fresh Spinach cache file and ignores the ordinary file that does not match the `spinach_*` cache pattern.
- Verifies `wipe_cache()` removes the Spinach cache file but keeps the unrelated ordinary file.
- Creates `spinach_old_dir` with a saved `payload.mat`, sets `tols.cache_mem=0`, and verifies `cacheman()` removes the matching out-of-date Spinach cache directory.
- Locates the shipped kernel cache directory as `fileparts(which('cacheman'))` and verifies the shipped cache records `st_product_table_2.mat` and `bos_product_table_2.mat` exist.
- Loads small cache-table records and checks dimensions: `st_product_table(2)`, `ist_product_table(2)`, and `bos_product_table(2)` each must return left and right tables of size `[4 4 4]`.
- Exercises the on-demand SLE operator cache: calls `sle_operators(1,2)`, and if the cache file `sle_operators_rank_1_int_2.mat` did not exist beforehand but was created by the call, deletes it. Verifies `Lx`, `Ly`, `Lz` are each `[10 10]`, `D` is a cell of size `[1 2]` with `D{1}` empty and `D{2}` of size `[5 5]`, and `space_basis` is `[10 3]`, matching the ten-function Wigner basis.
- Smoke-tests the read-only integrity sniffer: changes into `fileparts(which('sniff'))` with an `onCleanup` restore of the original directory, runs `sniff('none')` under `evalc`, and checks that it completes and prints either an all-clear or a fishy-file report.\n
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
