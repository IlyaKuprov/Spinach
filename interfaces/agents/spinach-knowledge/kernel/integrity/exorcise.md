# kernel/integrity/exorcise.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/integrity/exorcise.m`
- Signature: `exorcise(mode)`
- Total lines: 476

## Purpose

Searches Spinach distribution folders for any functions that do not conform to the house style. Opens the first one and complains to the console. Syntax: exorcise(mode)

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `is_block_start_token()`, `control_tokens()`, `strip_strings_and_comments()`, `is_block_end_token()`, `find_block_end()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(mode)`.
- Lines 28-29: List top level directories; implemented by `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 31-32: External package directories; implemented by `foreign_dirs={'jsonlab-1.5'}`.
- Lines 34-35: Get the directory trees; implemented by `mfiles=cell(numel(top_level,1))`.
- Lines 44-45: Process the files; implemented by `for k=1:numel(mfiles)`.
- Lines 47-48: Update the user; implemented by `disp(['inspecting file ' num2str(k) ' out of ' num2str(numel(mfiles)) ' '])`.
- Lines 50-51: Exclude foreign package directories; implemented by `if ~any(cellfun(@(x)contains(mfiles(k).folder,x),foreign_dirs))`.
- Lines 53-54: Read the file; implemented by `file_name=[mfiles(k).folder filesep mfiles(k).name]`.
- Lines 60-62: Identify top-level construct; implemented by `code_line_idx=find((~cellfun(@isempty,content))& (~startsWith(strtrim(content),'%')),1,'first')`.
- Lines 78-79: Apply exceptions; implemented by `grumbler_ex=false; header_ex=false`.
- Lines 96-97: Check for excessive white space; implemented by `blank_lines=cellfun(@isempty,content)`.
- Lines 104-105: Check for hard tab characters; implemented by `for m=1:numel(raw_content)`.
- Lines 111-112: Check for two line breaks at end of file; implemented by `file_text=fileread(file_name)`.
- Lines 117-118: Check that the grumbler is present; implemented by `if ~grumbler_ex`.
- Lines 133-135: Check that grumbler call is present in main body; implemented by `if (~grumbler_ex)&&is_function&&(~is_classdef)&&has_input_args&& (~contains(mfiles(k).folder,{'legacy','includes','overloads'}))`.
- Lines 137-138: Locate local function declarations; implemented by `function_lines=find(startsWith(strtrim(content),'function'))`.
- Lines 140-141: Find end of main function; implemented by `if numel(function_lines)>1`.
- Lines 147-148: Find grumble call in main function body; implemented by `grumble_call=false`.

### Control flow inferred from the code

- Line 37: `for` loop over `n=1:numel(top_level)`.
- Line 45: `for` loop over `k=1:numel(mfiles)`.
- Line 51: conditional branch on `~any(cellfun(@(x)contains(mfiles(k).folder,x),foreign_dirs))`.
- Line 63: conditional branch on `isempty(code_line_idx)`.
- Line 81: `for` loop over `m=1:numel(content)`.
- Line 82: conditional branch on `contains(content{m},'#NGRUM')`.
- Line 85: conditional branch on `contains(content{m},'#NHEAD')`.
- Line 88: conditional branch on `contains(content{m},'#NWIKI')`.
- Line 91: conditional branch on `contains(content{m},'#NORMOK')`.
- Line 98: `for` loop over `m=2:numel(blank_lines)`.
- Line 99: conditional branch on `blank_lines(m-1)&&blank_lines(m)`.
- Line 105: `for` loop over `m=1:numel(raw_content)`.
- Line 106: conditional branch on `contains(raw_content{m},sprintf('\t'))`.
- Line 113: conditional branch on `isempty(regexp(file_text,'(\r\n|\n){2}$','once'))`.

### Key state/data transformations

- Lines 29: computes `top_level` using `top_level={'kernel','interfaces','experiments','etc'}`.
- Lines 32: computes `foreign_dirs` using `foreign_dirs={'jsonlab-1.5'}`.
- Lines 35: computes `mfiles` using `mfiles=cell(numel(top_level,1))`.
- Lines 36: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-8))`.
- Lines 38-39: computes `mfiles{n}` using `mfiles{n}=dir([P '..' filesep '..' filesep top_level{n} filesep '**' filesep '*.m'])`.
- Lines 54: computes `file_name` using `file_name=[mfiles(k).folder filesep mfiles(k).name]`.
- Lines 55: computes `fid` using `fid=fopen(file_name,'r')`.
- Lines 56: computes `raw_content` using `raw_content=textscan(fid,'%s','Delimiter','\n','Whitespace','')`.
- Lines 57: computes `fclose(fid); raw_content` using `fclose(fid); raw_content=raw_content{1}`.
- Lines 58: computes `content` using `content=deblank(raw_content)`.
- Lines 61-62: computes `code_line_idx` using `code_line_idx=find((~cellfun(@isempty,content))& (~startsWith(strtrim(content),'%')),1,'first')`.
- Lines 64: computes `top_line` using `top_line=''`.
- Lines 65: computes `is_function` using `is_function=false`.
- Lines 66: computes `is_classdef` using `is_classdef=false`.
- Lines 67: computes `has_input_args` using `has_input_args=false`.
- Lines 68: computes `has_output_args` using `has_output_args=false`.
- Lines 73: computes `arg_block` using `arg_block=regexp(top_line,'\(([^)]*)\)','tokens','once')`.
- Lines 79: computes `grumbler_ex` using `grumbler_ex=false; header_ex=false`.

### Local helper functions

- Line 375: `is_block_start_token()` — `function answer=is_block_start_token(token)`. Extract control-flow tokens from a line
  - Representative operation: `answer=ismember(token,{'if','for','parfor','while', 'switch','try','spmd','function','classdef'})`.
  - Representative operation: `'switch','try','spmd','function','classdef'})`.
- Line 381: `control_tokens()` — `function [tokens,token_start,token_end,code_line]=control_tokens(line_text)`. Remove string literals and comments from a line
  - Representative operation: `code_line=lower(strip_strings_and_comments(line_text))`.
  - Representative operation: `[token_start,token_end,tokens]=regexp(code_line, '(?<![A-Za-z0-9_])(if|for|parfor|while|switch|try|spmd|function|classdef|otherwise|end)(?![A-Za-z0-9_])', 'start','end',…`.
- Line 389: `strip_strings_and_comments()` — `function code_line=strip_strings_and_comments(line_text)`.
  - Representative operation: `code_chars=[]; in_string=false; n=1`.
  - Representative operation: `while n<=numel(line_text)`.
- Line 419: `is_block_end_token()` — `function answer=is_block_end_token(code_line,start_idx,end_idx)`.
  - Representative operation: `prev_idx=find(~isspace(code_line(1:(start_idx-1))),1,'last')`.
  - Representative operation: `next_rel=find(~isspace(code_line((end_idx+1):end)),1,'first')`.
- Line 437: `find_block_end()` — `function end_line=find_block_end(content,start_line)`. Tokenise current line
  - Representative operation: `block_depth=1; end_line=[]`.
  - Representative operation: `for n=(start_line+1):numel(content)`.
- Line 462: `grumble()` — `function grumble(mode)`. It is my ambition to say in ten sentences what everyone else says in a whole book -what everyone else does not
  - Representative operation: `if ~ischar(mode)`.
  - Representative operation: `error('mode must be a character string.')`.

## Parameters / inputs

- mode -'online' checks the documentation
- Wiki for the corresponding page; 'offline'
- skips the Wiki check
- Any user contribution that this function has something to
- say about will either be brought under the house style, or
- rejected back to the user, depending on the amount of work
- involved. Always run this function before a commit if you
- have write access to Spinach repository.

## Implementation structure

- Searches Spinach distribution folders for any functions that
- do not conform to the house style. Opens the first one and
- complains to the console. Syntax:
- exorcise(mode)
- mode -'online' checks the documentation
- Wiki for the corresponding page; 'offline'
- skips the Wiki check
- Any user contribution that this function has something to
- say about will either be brought under the house style, or
- rejected back to the user, depending on the amount of work
- involved. Always run this function before a commit if you
- have write access to Spinach repository.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mfilename()`, `dir()`, `cell2mat()`, `mfiles()`, `randperm()`, `num2str()`, `any()`, `cellfun()`, `contains()`, `fopen()`, `textscan()`, `fclose()`, `deblank()`, `startsWith()`, `strtrim()`.
