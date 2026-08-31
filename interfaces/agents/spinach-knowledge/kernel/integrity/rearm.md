# kernel/integrity/rearm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/integrity/rearm.m`
- Signature: `rearm()`
- Total lines: 68

## Purpose

Rearms the sniffer database. The sniffer checks Spinach distribution .m files for any modifications that the user did since downloading Spinach. The function prints the list of files that have changed in any way since the in- ternal database has been rearmed. The purpose is to catch local modifications that the user may have made and forgotten about, that are causing some un- intended consequences elsewhere in Spinac

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: List top level directories; implemented by `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 20-21: List exceptions; implemented by `exceptions={}`.
- Lines 23-24: Get the directory trees; implemented by `mfiles=cell(numel(top_level,1))`.
- Lines 31-32: Get the table going; implemented by `smells={}`.
- Lines 34-35: Process the files; implemented by `for n=1:numel(top_level)`.
- Lines 39-40: Read the file; implemented by `file_name=[mfiles{n}(k).folder filesep mfiles{n}(k).name]`.
- Lines 45-46: Get the smell; implemented by `smells{end+1}=md5_hash([mfiles{n}(k).name md5_hash(content)])`.
- Lines 52-53: Overwrite the file; implemented by `delete([P 'smells.mat']); drawnow`.
- Lines 56-57: Inform the user; implemented by `disp('rearm: sniffer rearmed.')`.

### Control flow inferred from the code

- Line 26: `for` loop over `n=1:numel(top_level)`.
- Line 35: `for` loop over `n=1:numel(top_level)`.
- Line 36: `for` loop over `k=1:numel(mfiles{n})`.
- Line 37: conditional branch on `~ismember(mfiles{n}(k).name,exceptions)`.

### Key state/data transformations

- Lines 18: computes `top_level` using `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 21: computes `exceptions` using `exceptions={}`.
- Lines 24: computes `mfiles` using `mfiles=cell(numel(top_level,1))`.
- Lines 25: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-5))`.
- Lines 27-28: computes `mfiles{n}` using `mfiles{n}=dir([P '..' filesep '..' filesep top_level{n} filesep '**' filesep '*.m'])`.
- Lines 32: computes `smells` using `smells={}`.
- Lines 40: computes `file_name` using `file_name=[mfiles{n}(k).folder filesep mfiles{n}(k).name]`.
- Lines 41: computes `fid` using `fid=fopen(file_name,'r')`.
- Lines 42: computes `content` using `content=textscan(fid,'%s','Delimiter', '\n')`.
- Lines 43: computes `fclose(fid); content` using `fclose(fid); content=content{1}`.
- Lines 46: computes `smells{end+1}` using `smells{end+1}=md5_hash([mfiles{n}(k).name md5_hash(content)])`.

## Implementation structure

- Rearms the sniffer database. The sniffer checks Spinach
- distribution .m files for any modifications that the user
- did since downloading Spinach. The function prints the
- list of files that have changed in any way since the in-
- ternal database has been rearmed.
- The purpose is to catch local modifications that the user
- may have made and forgotten about, that are causing some un-
- intended consequences elsewhere in Spinach.
- List top level directories
- List exceptions
- Get the directory trees
- Get the table going

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mfilename()`, `dir()`, `ismember()`, `fopen()`, `textscan()`, `fclose()`, `md5_hash()`, `delete()`, `save()`.
