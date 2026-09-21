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
