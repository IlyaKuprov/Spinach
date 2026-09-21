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
