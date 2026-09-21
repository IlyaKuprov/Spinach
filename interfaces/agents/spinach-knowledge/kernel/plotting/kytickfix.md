# kernel/plotting/kytickfix.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/kytickfix.m`
- Signature: `kytickfix()`
- Total lines: 114

## Purpose

Switches Y axis tick labels to engineering notation by setting the numeric-ruler exponent to a multiple of three and installs a pan/zoom auto-updater. Syntax: kytickfix()

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `getappdata()`, `local_update()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- none

## Outputs

- updates the current axis system and installs an auto-updater

## Implementation structure

- Switches Y axis tick labels to engineering notation by setting
- the numeric-ruler exponent to a multiple of three and installs
- a pan/zoom auto-updater. Syntax:
- kytickfix()
- none
- updates the current axis system and installs an auto-updater
- Current Y axis ruler
- Check consistency
- Preserve any existing limit-change callback
- Install the automatic exponent updater
- Set the current exponent
- Local update callback

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isappdata()`, `getappdata()`, `isfield()`, `isequal()`, `local_callback()`, `setappdata()`, `local_update()`, `ishandle()`, `iscell()`, `feval()`, `ischar()`, `isstring()`, `evalin()`, `char()`, `strcmp()`.
