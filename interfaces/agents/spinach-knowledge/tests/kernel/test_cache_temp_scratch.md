# tests/kernel/test_cache_temp_scratch.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_cache_temp_scratch.m`
- Signature: `result=test_cache_temp_scratch()`
- Total lines: 120

## Purpose

Tests cache management in a temporary scratch directory. Syntax: result=test_cache_temp_scratch()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_remove_dir()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks cacheman() and wipe_cache() against a temporary scratch
- directory, loads shipped small cache tables without generating new kernel
- cache files, and smoke-tests the read-only sniff() integrity pass.

## Implementation structure

- Tests cache management in a temporary scratch directory. Syntax:
- result=test_cache_temp_scratch()
- result -regression test result with explanatory messages
- The test checks cacheman() and wipe_cache() against a temporary scratch
- directory, loads shipped small cache tables without generating new kernel
- cache files, and smoke-tests the read-only sniff() integrity pass.
- Announce the test target
- State the cache-management target of the test
- Create an isolated scratch directory and arrange cleanup
- Build a minimal Spinach object pointing at the temporary scratch directory
- Ensure cacheman() uses a small process pool instead of auto-starting a large one
- Create Spinach and non-Spinach scratch files

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `cacheman()`, `wipe_cache()`, `tempname()`, `mkdir()`, `onCleanup()`, `local_remove_dir()`, `gcp()`, `parpool()`, `fullfile()`, `save()`, `test_true()`, `exist()`, `fileparts()`, `which()`, `st_product_table()`.
