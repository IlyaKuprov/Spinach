# tests/run_tests.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/run_tests.m`
- Signature: `results=run_tests(varargin)`
- Total lines: 94

## Purpose

Runs the Spinach regression test suite. Syntax: results=run_tests(varargin)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Add the test library to the path; implemented by `root_dir=fileparts(mfilename('fullpath'))`.
- Lines 25-26: Add the Spinach production directories to the path; implemented by `spinach_root=fileparts(root_dir)`.
- Lines 32-33: Parse options; implemented by `options=test_options(varargin{:})`.
- Lines 35-36: Get the manifest; implemented by `manifest=test_manifest()`.
- Lines 38-39: Apply substring filter; implemented by `if ~isempty(options.pattern)`.
- Lines 45-47: Preallocate result array; implemented by `results=struct('id',{},'name',{},'purpose',{},'status',{},'elapsed',{}, 'messages',{},'error',{})`.
- Lines 49-50: Run the tests; implemented by `for n=1:numel(manifest)`.
- Lines 81-82: Summarise outcomes; implemented by `n_pass=nnz(strcmp({results.status},'PASS'))`.

### Control flow inferred from the code

- Line 39: conditional branch on `~isempty(options.pattern)`.
- Line 50: `for` loop over `n=1:numel(manifest)`.
- Line 67: conditional branch on `options.verbose`.
- Line 69: `for` loop over `k=1:numel(result.messages)`.
- Line 72: conditional branch on `strcmp(result.status,'FAIL')`.
- Line 76: conditional branch on `options.stop_on_fail&&strcmp(result.status,'FAIL')`.
- Line 85: conditional branch on `n_fail>0`.
- Line 87: `for` loop over `n=1:numel(failed)`.

### Key state/data transformations

- Lines 19: computes `root_dir` using `root_dir=fileparts(mfilename('fullpath'))`.
- Lines 26: computes `spinach_root` using `spinach_root=fileparts(root_dir)`.
- Lines 33: computes `options` using `options=test_options(varargin{:})`.
- Lines 36: computes `manifest` using `manifest=test_manifest()`.
- Lines 40-41: computes `keep` using `keep=contains({manifest.id},options.pattern)| contains({manifest.name},options.pattern)`.
- Lines 46-47: computes `results` using `results=struct('id',{},'name',{},'purpose',{},'status',{},'elapsed',{}, 'messages',{},'error',{})`.
- Lines 53: computes `result` using `result=feval(manifest(n).function)`.
- Lines 54: computes `result.status` using `result.status='PASS'`.
- Lines 55: computes `result.elapsed` using `result.elapsed=toc`.
- Lines 56: computes `result.error` using `result.error=''`.
- Lines 58: computes `result.id` using `result.id=manifest(n).id`.
- Lines 59: computes `result.name` using `result.name=manifest(n).name`.
- Lines 60: computes `result.purpose` using `result.purpose=''`.
- Lines 63: computes `result.messages` using `result.messages={}`.
- Lines 66: computes `results(end+1)` using `results(end+1)=result`.
- Lines 82: computes `n_pass` using `n_pass=nnz(strcmp({results.status},'PASS'))`.
- Lines 83: computes `n_fail` using `n_fail=nnz(strcmp({results.status},'FAIL'))`.
- Lines 86: computes `failed` using `failed=results(strcmp({results.status},'FAIL'))`.

## Parameters / inputs

- varargin -name-value options: 'pattern', 'verbose', and
- 'stop_on_fail'

## Outputs

- results -structure array with test outcomes and messages

## Implementation structure

- Runs the Spinach regression test suite. Syntax:
- results=run_tests(varargin)
- varargin -name-value options: 'pattern', 'verbose', and
- 'stop_on_fail'
- results -structure array with test outcomes and messages
- Add the test library to the path
- Add the Spinach production directories to the path
- Parse options
- Get the manifest
- Apply substring filter
- Preallocate result array
- Run the tests

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fileparts()`, `mfilename()`, `addpath()`, `fullfile()`, `genpath()`, `test_options()`, `test_manifest()`, `contains()`, `manifest()`, `feval()`, `results()`, `strcmp()`, `nnz()`, `failed()`.
