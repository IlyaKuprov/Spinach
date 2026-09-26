# kernel/cache/wipe_cache.m

- Signature: `wipe_cache(spin_system)`

## Purpose

Forces a wipe of the Spinach cache folder. Syntax: wipe_cache(spin_system)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach object with information (stored
- in spin_system.sys.scratch) about the
- cache folder location, use bootstrap()
- to get the default object
- Output:
- an attempt is made to delete all all Spinach-specific
- files in spin_system.sys.scratch; this would fail qui-
- etly if file system permissions are insufficient

## Implementation structure

- Forces a wipe of the Spinach cache folder. Syntax:
- wipe_cache(spin_system)
- spin_system -Spinach object with information (stored
- in spin_system.sys.scratch) about the
- cache folder location, use bootstrap()
- to get the default object
- Output:
- an attempt is made to delete all all Spinach-specific
- files in spin_system.sys.scratch; this would fail qui-
- etly if file system permissions are insufficient
- Defaults for command line calls
- Check consistency
