# kernel/integrity/existentials.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/integrity/existentials.m`
- Signature: `existentials()`
- Total lines: 130

## Purpose

Kernel integrity control. Checks for collisions between Spinach functions and anything else that the user may have installed or written in the current Matlab instance. Also checks for any fi- les that are not visible to Matlab because the corresponding di- rectory is not on the path. Collisions of function names and path problems are the most fre- quent support topic at the forum.

## Physical / mathematical content

- Integrity-control utilities. These files check distribution state, path collisions, style conformance, sniffer databases, and other safeguards that protect Spinach reproducibility.

## Numerical / algorithmic content

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
