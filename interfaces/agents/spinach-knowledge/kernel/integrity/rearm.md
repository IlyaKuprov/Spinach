# kernel/integrity/rearm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/integrity/rearm.m`
- Signature: `rearm()`
- Total lines: 68

## Purpose

Rearms the sniffer database. The sniffer checks Spinach distribution .m files for any modifications that the user did since downloading Spinach. The function prints the list of files that have changed in any way since the in- ternal database has been rearmed. The purpose is to catch local modifications that the user may have made and forgotten about, that are causing some un- intended consequences elsewhere in Spinac

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

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
