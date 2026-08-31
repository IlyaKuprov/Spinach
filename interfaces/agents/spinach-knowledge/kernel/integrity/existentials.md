# kernel/integrity/existentials.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/integrity/existentials.m`
- Signature: `existentials()`
- Total lines: 130

## Purpose

Kernel integrity control. Checks for collisions between Spinach functions and anything else that the user may have installed or written in the current Matlab instance. Also checks for any fi- les that are not visible to Matlab because the corresponding di- rectory is not on the path. Collisions of function names and path problems are the most fre- quent support topic at the forum.

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Do not run inside parallel pools; implemented by `if isworkernode, return; end`.
- Lines 20-21: Inform the user; implemented by `disp('Running startup checks ')`.
- Lines 23-29: ########################################## NO, IT WILL NOT MAGICALLY START WORKING % IF YOU COMMENT ANY OF THIS OUT % ##########################################; implemented by `if isMATLABReleaseOlderThan('R2024b')`.
- Lines 28-29: Global Matlab release; implemented by `if isMATLABReleaseOlderThan('R2024b')`.
- Lines 33-34: Existential toolboxes; implemented by `if ~exist([matlabroot filesep 'toolbox' filesep 'parallel'],'dir')`.
- Lines 56-57: List top level directories; implemented by `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 59-60: Get the directory trees; implemented by `mfiles=cell(numel(top_level,1))`.
- Lines 67-68: Process the files; implemented by `for n=1:numel(top_level)`.
- Lines 71-72: Get the full name as per Spinach; implemented by `file_name_spinach=[mfiles{n}(k).folder filesep mfiles{n}(k).name]`.
- Lines 74-75: Get the full name as per Matlab; implemented by `file_name_matlab=which(mfiles{n}(k).name)`.
- Lines 77-80: Check for collisions, except in overloads; implemented by `if (~isempty(file_name_matlab))&& (~strcmp(file_name_spinach,file_name_matlab))&& (~contains(mfiles{n}(k).folder,'overloads'))`.
- Lines 91-92: Check for missing files; implemented by `if isempty(file_name_matlab)`.
- Lines 106-107: Windows; implemented by `if ispc`.
- Lines 109-110: Find out which file system Spinach volume has; implemented by `own_disk=mfilename('fullpath'); own_disk=own_disk(1:3)`.
- Lines 113-114: Warn about non-NTFS volumes; implemented by `if ~strcmp(char(own_disk.DriveFormat),'NTFS')`.

### Control flow inferred from the code

- Line 18: conditional branch on `isworkernode, return; end`.
- Line 29: conditional branch on `isMATLABReleaseOlderThan('R2024b')`.
- Line 34: conditional branch on `~exist([matlabroot filesep 'toolbox' filesep 'parallel'],'dir')`.
- Line 37: conditional branch on `~exist([matlabroot filesep 'toolbox' filesep 'nnet'],'dir')`.
- Line 40: conditional branch on `~exist([matlabroot filesep 'toolbox' filesep 'rl'],'dir')`.
- Line 43: conditional branch on `~exist([matlabroot filesep 'toolbox' filesep 'optim'],'dir')`.
- Line 46: conditional branch on `~exist([matlabroot filesep 'toolbox' filesep 'stats'],'dir')`.
- Line 49: conditional branch on `~exist([matlabroot filesep 'toolbox' filesep 'map'],'dir')`.
- Line 52: conditional branch on `~exist([matlabroot filesep 'toolbox' filesep 'aero'],'dir')`.
- Line 62: `for` loop over `n=1:numel(top_level)`.
- Line 68: `for` loop over `n=1:numel(top_level)`.
- Line 69: `for` loop over `k=1:numel(mfiles{n})`.
- Line 78: conditional branch on `(~isempty(file_name_matlab))&&`.
- Line 92: conditional branch on `isempty(file_name_matlab)`.

### Key state/data transformations

- Lines 57: computes `top_level` using `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 60: computes `mfiles` using `mfiles=cell(numel(top_level,1))`.
- Lines 61: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-13))`.
- Lines 63-64: computes `mfiles{n}` using `mfiles{n}=dir([P filesep '..' filesep '..' filesep top_level{n} filesep '**' filesep '*.m'])`.
- Lines 72: computes `file_name_spinach` using `file_name_spinach=[mfiles{n}(k).folder filesep mfiles{n}(k).name]`.
- Lines 75: computes `file_name_matlab` using `file_name_matlab=which(mfiles{n}(k).name)`.
- Lines 110: computes `own_disk` using `own_disk=mfilename('fullpath'); own_disk=own_disk(1:3)`.

## Implementation structure

- Kernel integrity control. Checks for collisions between Spinach
- functions and anything else that the user may have installed or
- written in the current Matlab instance. Also checks for any fi-
- les that are not visible to Matlab because the corresponding di-
- rectory is not on the path.
- Collisions of function names and path problems are the most fre-
- quent support topic at the forum.
- Do not run inside parallel pools
- Inform the user
- ##########################################
- NO, IT WILL NOT MAGICALLY START WORKING %
- IF YOU COMMENT ANY OF THIS OUT %

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isMATLABReleaseOlderThan()`, `exist()`, `mfilename()`, `dir()`, `which()`, `strcmp()`, `contains()`, `own_disk()`, `char()`.
