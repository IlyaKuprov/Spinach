# kernel/integrity/sniff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/integrity/sniff.m`
- Signature: `sniff(action)`
- Total lines: 118

## Purpose

Kernel integrity control. Checks Spinach distribution .m files for any modifications that the user did since downloading Spi- nach. The function prints the list of files that have changed in any way since the internal database has been rearmed. The purpose is to catch local modifications that the user may have made and forgotten about, that are causing some unintend- ed consequences elsewhere in Spinach. Syntax: snif

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Default is to take no action; implemented by `if ~exist('action','var'), action='none'; end`.
- Lines 26-27: Check consistency; implemented by `grumble(action)`.
- Lines 29-30: List top level directories; implemented by `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 32-33: List exceptions; implemented by `exceptions={}`.
- Lines 35-36: Get the directory trees; implemented by `mfiles=cell(numel(top_level,1))`.
- Lines 43-44: Read the smells table; implemented by `load('smells.mat','smells'); all_good=true()`.
- Lines 46-47: Start line counters; implemented by `code_line_count=0`.
- Lines 50-51: Process the files; implemented by `for n=1:numel(top_level)`.
- Lines 55-56: Read the file; implemented by `file_name=[mfiles{n}(k).folder filesep mfiles{n}(k).name]`.
- Lines 61-62: Get the smell; implemented by `smell=md5_hash([mfiles{n}(k).name md5_hash(content)])`.
- Lines 64-65: Drop blank lines; implemented by `content(cellfun(@(x)isempty(deblank(x)),content))=[]`.
- Lines 67-69: Count comment lines; implemented by `comm_line_count=comm_line_count+ nnz(cellfun(@(x)strcmp('%',x(1)),content))`.
- Lines 71-72: Drop comment lines; implemented by `content(cellfun(@(x)strcmp('%',x(1)),content))=[]`.
- Lines 74-75: Count code lines; implemented by `code_line_count=code_line_count+numel(content)`.
- Lines 77-78: Check for a match; implemented by `if ~ismember(smell,smells)`.
- Lines 80-81: Previously unseen file; implemented by `disp(['smells fishy: ' file_name]); all_good=false()`.
- Lines 83-84: Open the file; implemented by `if strcmp(action,'open'), edit(file_name); end`.
- Lines 92-93: Give an all-clear if appropriate; implemented by `if all_good`.

### Control flow inferred from the code

- Line 24: conditional branch on `~exist('action','var'), action='none'; end`.
- Line 38: `for` loop over `n=1:numel(top_level)`.
- Line 51: `for` loop over `n=1:numel(top_level)`.
- Line 52: `for` loop over `k=1:numel(mfiles{n})`.
- Line 53: conditional branch on `~ismember(mfiles{n}(k).name,exceptions)`.
- Line 78: conditional branch on `~ismember(smell,smells)`.
- Line 84: conditional branch on `strcmp(action,'open'), edit(file_name); end`.
- Line 93: conditional branch on `all_good`.

### Key state/data transformations

- Lines 30: computes `top_level` using `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 33: computes `exceptions` using `exceptions={}`.
- Lines 36: computes `mfiles` using `mfiles=cell(numel(top_level,1))`.
- Lines 37: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-5))`.
- Lines 39-40: computes `mfiles{n}` using `mfiles{n}=dir([P filesep '..' filesep '..' filesep top_level{n} filesep '**' filesep '*.m'])`.
- Lines 44: computes `load('smells.mat','smells'); all_good` using `load('smells.mat','smells'); all_good=true()`.
- Lines 47: computes `code_line_count` using `code_line_count=0`.
- Lines 48: computes `comm_line_count` using `comm_line_count=0`.
- Lines 56: computes `file_name` using `file_name=[mfiles{n}(k).folder filesep mfiles{n}(k).name]`.
- Lines 57: computes `fid` using `fid=fopen(file_name,'r')`.
- Lines 58: computes `content` using `content=textscan(fid,'%s','Delimiter', '\n')`.
- Lines 59: computes `fclose(fid); content` using `fclose(fid); content=content{1}`.
- Lines 62: computes `smell` using `smell=md5_hash([mfiles{n}(k).name md5_hash(content)])`.
- Lines 65: computes `content(cellfun(@(x)isempty(deblank(x)),content))` using `content(cellfun(@(x)isempty(deblank(x)),content))=[]`.
- Lines 72: computes `content(cellfun(@(x)strcmp('%',x(1)),content))` using `content(cellfun(@(x)strcmp('%',x(1)),content))=[]`.
- Lines 81: computes `disp(['smells fishy: ' file_name]); all_good` using `disp(['smells fishy: ' file_name]); all_good=false()`.

### Local helper functions

- Line 102: `grumble()` — `function grumble(action)`. Glen had an unusual introduction to the Arctic in 1932. He thought he had accepted a friend's invitation to a debutante
  - Representative operation: `if (~ischar(action))||(~ismember(action,{'none','open'}))`.
  - Representative operation: `error('action must be ''none'' or ''open''.')`.

## Parameters / inputs

- action -'none' prints the names of fishy files to
- the console, 'open' opens them

## Implementation structure

- Kernel integrity control. Checks Spinach distribution .m files
- for any modifications that the user did since downloading Spi-
- nach. The function prints the list of files that have changed
- in any way since the internal database has been rearmed.
- The purpose is to catch local modifications that the user may
- have made and forgotten about, that are causing some unintend-
- ed consequences elsewhere in Spinach. Syntax:
- sniff(action)
- action -'none' prints the names of fishy files to
- the console, 'open' opens them
- Default is to take no action
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `mfilename()`, `dir()`, `load()`, `true()`, `ismember()`, `fopen()`, `textscan()`, `fclose()`, `md5_hash()`, `content()`, `cellfun()`, `deblank()`, `nnz()`, `strcmp()`.
