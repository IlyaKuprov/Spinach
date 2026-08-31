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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Temporary scratch cache management\n')`.
- Lines 20-23: State the cache-management target of the test; implemented by `result=new_test_result('kernel/cache_temp_scratch', 'Temporary scratch cache management', 'cacheman() and wipe_cache() must only affect Spinach cache records in the conf…`.
- Lines 25-26: Create an isolated scratch directory and arrange cleanup; implemented by `scratch=tempname(tempdir)`.
- Lines 30-31: Build a minimal Spinach object pointing at the temporary scratch directory; implemented by `spin_system.sys.output='hush'`.
- Lines 37-38: Ensure cacheman() uses a small process pool instead of auto-starting a large one; implemented by `current_pool=gcp('nocreate')`.
- Lines 43-44: Create Spinach and non-Spinach scratch files; implemented by `payload=1`.
- Lines 50-51: Check that a fresh Spinach cache file is kept by a long retention horizon; implemented by `cacheman(spin_system)`.
- Lines 57-58: Check that wipe_cache() deletes Spinach cache files but not unrelated files; implemented by `wipe_cache(spin_system)`.
- Lines 64-65: Check that cacheman() removes matching cache directories under a zero horizon; implemented by `spinach_dir=fullfile(scratch,'spinach_old_dir')`.
- Lines 73-74: Locate the shipped kernel cache directory; implemented by `cache_root=fileparts(which('cacheman'))`.
- Lines 76-78: Verify that small shipped cache records exist before loading them; implemented by `result=test_true(result,'st_product_table shipped cache',exist(fullfile(cache_root,'st_product_table_2.mat'),'file')==2, 'the shipped two-level ST product table cache sh…`.
- Lines 84-85: Load small cache-table records and check their dimensions; implemented by `[st_left,st_right]=st_product_table(2)`.
- Lines 100-101: Smoke-test the read-only integrity sniffer without requiring a pristine tree; implemented by `old_dir=pwd`.

### Control flow inferred from the code

- Line 39: conditional branch on `isempty(current_pool)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/cache_temp_scratch', 'Temporary scratch cache management', 'cacheman() and wipe_cache() must only affect Spinach cache records in the conf…`.
- Lines 26: computes `scratch` using `scratch=tempname(tempdir)`.
- Lines 28: computes `cleanup_obj` using `cleanup_obj=onCleanup(@()local_remove_dir(scratch))`.
- Lines 31: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 32: computes `spin_system.sys.scratch` using `spin_system.sys.scratch=scratch`.
- Lines 33: computes `spin_system.sys.enable` using `spin_system.sys.enable={}`.
- Lines 34: computes `spin_system.sys.disable` using `spin_system.sys.disable={}`.
- Lines 35: computes `spin_system.tols.cache_mem` using `spin_system.tols.cache_mem=365`.
- Lines 38: computes `current_pool` using `current_pool=gcp('nocreate')`.
- Lines 44: computes `payload` using `payload=1`.
- Lines 45: computes `spinach_file` using `spinach_file=fullfile(scratch,'spinach_fresh.mat')`.
- Lines 46: computes `ordinary_file` using `ordinary_file=fullfile(scratch,'ordinary_file.mat')`.
- Lines 65: computes `spinach_dir` using `spinach_dir=fullfile(scratch,'spinach_old_dir')`.
- Lines 74: computes `cache_root` using `cache_root=fileparts(which('cacheman'))`.
- Lines 85: computes `[st_left,st_right]` using `[st_left,st_right]=st_product_table(2)`.
- Lines 88: computes `[ist_left,ist_right]` using `[ist_left,ist_right]=ist_product_table(2)`.
- Lines 91: computes `[bos_left,bos_right]` using `[bos_left,bos_right]=bos_product_table(2)`.
- Lines 94: computes `[Lx,Ly,Lz,D,space_basis]` using `[Lx,Ly,Lz,D,space_basis]=sle_operators(1,2)`.

### Local helper functions

- Line 111: `local_remove_dir()` — `function local_remove_dir(path_name)`. Remove a temporary directory if it still exists
  - Representative operation: `if exist(path_name,'dir')`.
  - Representative operation: `rmdir(path_name,'s')`.

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
