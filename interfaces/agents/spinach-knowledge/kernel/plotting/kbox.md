# kernel/plotting/kbox.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/kbox.m`
- Signature: `kbox() % #NGRUM`
- Total lines: 133

## Purpose

Creates a tickless boxed frame around the current axes using ordinary line objects that live in the same data space as the plot. This is needed for spectrograms be- cause Matlab has dumb plotting defaults. Syntax: kbox()

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `getappdata()`, `local_seg()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- creates or updates a tickless axis box
- in the current axes

## Implementation structure

- Creates a tickless boxed frame around the current axes
- using ordinary line objects that live in the same data
- space as the plot. This is needed for spectrograms be-
- cause Matlab has dumb plotting defaults. Syntax:
- kbox()
- creates or updates a tickless axis box
- in the current axes
- Get the current axis object
- Remove any previous Spinach box overlay
- Delete orphaned overlay axes from the rejected implementation
- Create the box line in the plot axes
- Exclude the box line from autoscaling

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isappdata()`, `getappdata()`, `rmappdata()`, `isfield()`, `cellfun()`, `ishandle()`, `delete()`, `findall()`, `get()`, `set()`, `line()`, `isprop()`, `findprop()`, `addlistener()`, `local_update()`, `setappdata()`.
