# tests/kernel/test_cache_temp_scratch.m

- Signature: `result=test_cache_temp_scratch()`

## Purpose

Tests cache management in a temporary scratch directory. Syntax: result=test_cache_temp_scratch()

## Physical / mathematical content

## Numerical / algorithmic content

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
