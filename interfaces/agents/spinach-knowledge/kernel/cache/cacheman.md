# kernel/cache/cacheman.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/cacheman.m`
- Signature: `cacheman(spin_system) %#NHEAD`
- Total lines: 103

## Purpose

Cache management heuristics. Looks after the scratch folder and prevents it from filling up the disk. Do not call directly. The function inspects the scratch folder and deletes any files that are older than the threshold (default is 365 days) speci- fied in spin_system.tols.cache_mem field.

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Check consistency; implemented by `grumble(spin_system)`.
- Lines 17-18: Get parallel pool directory; implemented by `current_pool=gcp('nocreate')`.
- Lines 25-26: Calculate the time horizon; implemented by `time_horizon=now-spin_system.tols.cache_mem`.
- Lines 35-36: Look into the scratch directory; implemented by `dir_cont=dir([spin_system.sys.scratch filesep 'spinach_*'])`.
- Lines 41-42: Delete anything that is out of date; implemented by `n_files_gone=0; n_dirs_gone=0`.
- Lines 61-62: Report to the user; implemented by `if n_files_gone>0`.

### Control flow inferred from the code

- Line 19: conditional branch on `~isempty(current_pool)`.
- Line 37: conditional branch on `any([dir_cont.datenum]<time_horizon)`.
- Line 43: `parfor` loop over `n=1:numel(dir_cont)`.
- Line 44: conditional branch on `dir_cont(n).datenum<time_horizon`.
- Line 45: conditional branch on `dir_cont(n).isdir`.
- Line 48: conditional branch on `~strcmp(dir_name,pool_dir)`.
- Line 62: conditional branch on `n_files_gone>0`.
- Line 65: conditional branch on `n_dirs_gone>0`.

### Key state/data transformations

- Lines 18: computes `current_pool` using `current_pool=gcp('nocreate')`.
- Lines 20: computes `pool_dir` using `pool_dir=current_pool.Cluster.JobStorageLocation`.
- Lines 26: computes `time_horizon` using `time_horizon=now-spin_system.tols.cache_mem`.
- Lines 29: computes `test` using `test=1; save([spin_system.sys.scratch filesep 'test.mat'],'test')`.
- Lines 36: computes `dir_cont` using `dir_cont=dir([spin_system.sys.scratch filesep 'spinach_*'])`.
- Lines 42: computes `n_files_gone` using `n_files_gone=0; n_dirs_gone=0`.
- Lines 47: computes `dir_name` using `dir_name=[dir_cont(n).folder filesep dir_cont(n).name]`.
- Lines 49: computes `rmdir(dir_name,'s'); n_dirs_gone` using `rmdir(dir_name,'s'); n_dirs_gone=n_dirs_gone+1`.

### Local helper functions

- Line 72: `grumble()` — `function grumble(spin_system)`.
  - Representative operation: `if (~isfield(spin_system,'sys'))||(~isfield(spin_system.sys,'scratch'))`.
  - Representative operation: `error('the spin_system object does not specify scratch location.')`.

## Implementation structure

- Cache management heuristics. Looks after the scratch folder and
- prevents it from filling up the disk. Do not call directly.
- The function inspects the scratch folder and deletes any files
- that are older than the threshold (default is 365 days) speci-
- fied in spin_system.tols.cache_mem field.
- Check consistency
- Get parallel pool directory
- Calculate the time horizon
- Look into the scratch directory
- Delete anything that is out of date
- Report to the user
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gcp()`, `save()`, `delete()`, `dir()`, `any()`, `report()`, `dir_cont()`, `strcmp()`, `rmdir()`, `num2str()`, `isfield()`, `ischar()`, `isscalar()`, `exist()`.
