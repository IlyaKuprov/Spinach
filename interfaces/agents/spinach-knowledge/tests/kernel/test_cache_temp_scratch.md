# tests/kernel/test_cache_temp_scratch.m

- Signature: `result=test_cache_temp_scratch()`

## Purpose

Tests cache management routines in an isolated temporary scratch directory.

## Numerical / algorithmic content

With `tols.cache_mem=365`, the test checks that `cacheman()` keeps a fresh `spinach_*` file and ignores an ordinary file. It checks that `wipe_cache()` removes the Spinach file while preserving the ordinary file, and that `cacheman()` removes a matching cache directory when `tols.cache_mem=0`. It also checks shipped small cache records and product-table dimensions, generates SLE operators, and deletes the SLE cache record only if this test created it.

## Parameters / inputs

The test uses a temporary `sys.scratch` directory, `sys.output='hush'`, empty `sys.enable` and `sys.disable`, and a one-process `Processes` pool only if no pool exists.

## Outputs

- `result` — regression test result with explanatory messages.
- The test smoke-tests `sniff('none')` for nonempty output without requiring a pristine tree.
