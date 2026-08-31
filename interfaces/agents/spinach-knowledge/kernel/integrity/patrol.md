# kernel/integrity/patrol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/integrity/patrol.m`
- Signature: `patrol(test_subject)`
- Total lines: 125

## Purpose

This function runs contiouously on one of our servers, its purpose is to catch any unintended consequences before they propagate too far down the development chain. Syntax: patrol(test_subject)

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Set default and check consistency; implemented by `if ~exist('test_subject','var')`.
- Lines 30-31: List exceptions; implemented by `exceptions={}`.
- Lines 33-34: Shuffle the RNG; implemented by `rng('shuffle')`.
- Lines 36-37: Get the directory tree; implemented by `P=mfilename('fullpath'); P=P(1:(end-6))`.
- Lines 41-42: Find relevant files; implemented by `if ~isempty(test_subject)`.
- Lines 44-45: No file is initially relevant; implemented by `relevant_file_mask=false(numel(mfiles),1)`.
- Lines 47-48: Loop over the file list; implemented by `for n=1:numel(mfiles)`.
- Lines 50-51: Get full name; implemented by `full_name=[mfiles(n).folder filesep mfiles(n).name]`.
- Lines 53-54: Read the file; implemented by `fid=fopen(full_name,'r')`.
- Lines 58-59: Files that mention the subject are relevant; implemented by `for m=1:numel(content)`.
- Lines 70-71: Test all example files; implemented by `relevant_file_mask=true(numel(mfiles),1)`.
- Lines 75-76: Update the file list; implemented by `mfiles=mfiles(relevant_file_mask)`.
- Lines 78-79: Inform the user; implemented by `if isempty(mfiles)`.
- Lines 85-86: Enter the main loop; implemented by `hashes_match=true()`.
- Lines 89-90: Pick a random file; implemented by `n=randi(numel(mfiles))`.
- Lines 92-93: Build the file name; implemented by `file_name=[mfiles(n).folder filesep mfiles(n).name]`.
- Lines 96-97: Check that Matlab's syntax checker is on green; implemented by `if ~isempty(checkcode(file_name))`.
- Lines 101-102: Run the file; implemented by `if ~ismember(mfiles(n).name,exceptions)`.

### Control flow inferred from the code

- Line 25: conditional branch on `~exist('test_subject','var')`.
- Line 42: conditional branch on `~isempty(test_subject)`.
- Line 48: `for` loop over `n=1:numel(mfiles)`.
- Line 59: `for` loop over `m=1:numel(content)`.
- Line 60: conditional branch on `contains(content{m},test_subject)||`.
- Line 79: conditional branch on `isempty(mfiles)`.
- Line 87: `while` loop over `hashes_match`.
- Line 97: conditional branch on `~isempty(checkcode(file_name))`.
- Line 102: conditional branch on `~ismember(mfiles(n).name,exceptions)`.

### Key state/data transformations

- Lines 26: computes `test_subject` using `test_subject=''`.
- Lines 31: computes `exceptions` using `exceptions={}`.
- Lines 37: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-6))`.
- Lines 38-39: computes `mfiles` using `mfiles=dir([P '..' filesep '..' filesep 'examples' filesep '**' filesep '*.m'])`.
- Lines 45: computes `relevant_file_mask` using `relevant_file_mask=false(numel(mfiles),1)`.
- Lines 51: computes `full_name` using `full_name=[mfiles(n).folder filesep mfiles(n).name]`.
- Lines 54: computes `fid` using `fid=fopen(full_name,'r')`.
- Lines 55: computes `content` using `content=textscan(fid,'%s','Delimiter','\n')`.
- Lines 56: computes `fclose(fid); content` using `fclose(fid); content=content{1}`.
- Lines 62: computes `relevant_file_mask(n)` using `relevant_file_mask(n)=true()`.
- Lines 86: computes `hashes_match` using `hashes_match=true()`.
- Lines 90: computes `n` using `n=randi(numel(mfiles))`.
- Lines 93: computes `file_name` using `file_name=[mfiles(n).folder filesep mfiles(n).name]`.

### Local helper functions

- Line 116: `grumble()` — `function grumble(test_subject)`. No snowflake in an avalanche ever feels responsible. Stanislav Lec
  - Representative operation: `if ~ischar(test_subject)`.
  - Representative operation: `error('test_subject must be a character string')`.

## Parameters / inputs

- test_subject -a character string; if it occurs
- anywhere within the example file
- path, that file is included into
- the patrol run

## Outputs

- whatever the individual examples return

## Implementation structure

- This function runs contiouously on one of our servers, its
- purpose is to catch any unintended consequences before they
- propagate too far down the development chain. Syntax:
- patrol(test_subject)
- test_subject -a character string; if it occurs
- anywhere within the example file
- path, that file is included into
- the patrol run
- whatever the individual examples return
- Set default and check consistency
- List exceptions
- Shuffle the RNG

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `rng()`, `mfilename()`, `dir()`, `false()`, `mfiles()`, `fopen()`, `textscan()`, `fclose()`, `contains()`, `relevant_file_mask()`, `true()`, `num2str()`, `randi()`, `checkcode()`.
