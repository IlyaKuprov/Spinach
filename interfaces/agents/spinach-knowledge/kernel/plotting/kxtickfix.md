# kernel/plotting/kxtickfix.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/kxtickfix.m`
- Signature: `kxtickfix()`
- Total lines: 115

## Purpose

Switches X axis tick labels to engineering notation by setting the numeric-ruler exponent to a multiple of three and installs a pan/zoom auto-updater. Syntax: kxtickfix()

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `getappdata()`, `local_update()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Current X axis ruler; implemented by `ax=gca; ruler=ax.XAxis`.
- Lines 24-25: Check consistency; implemented by `grumble(ruler)`.
- Lines 27-28: Preserve any existing limit-change callback; implemented by `app_key='SpinachKXTickFix'`.
- Lines 40-41: Install the automatic exponent updater; implemented by `callback=@(source,event)local_callback(source,event,ax,app_key)`.
- Lines 46-47: Set the current exponent; implemented by `local_update(ruler)`.

### Control flow inferred from the code

- Line 29: conditional branch on `isappdata(ax,app_key)`.
- Line 31: conditional branch on `isfield(state,'callback')&&isequal(ruler.LimitsChangedFcn,state.callback)`.

### Key state/data transformations

- Lines 22: computes `ax` using `ax=gca; ruler=ax.XAxis`.
- Lines 28: computes `app_key` using `app_key='SpinachKXTickFix'`.
- Lines 30: computes `state` using `state=getappdata(ax,app_key)`.
- Lines 32: computes `old_fcn` using `old_fcn=state.old_fcn`.
- Lines 41: computes `callback` using `callback=@(source,event)local_callback(source,event,ax,app_key)`.
- Lines 42: computes `state.callback` using `state.callback=callback; state.old_fcn=old_fcn`.
- Lines 44: computes `ruler.LimitsChangedFcn` using `ruler.LimitsChangedFcn=callback`.

### Local helper functions

- Line 52: `local_callback()` — `function local_callback(ruler,event,ax,app_key)`. Ignore deleted graphics objects
  - Representative operation: `if ~ishandle(ruler)||~ishandle(ax)||~isappdata(ax,app_key), return, end`.
  - Representative operation: `local_update(ruler)`.
- Line 75: `local_update()` — `function local_update(ruler)`. Ignore non-linear axis state after installation
  - Representative operation: `if ~strcmp(ruler.Scale,'linear'), return, end`.
  - Representative operation: `ticks=ruler.TickValues`.
- Line 100: `grumble()` — `function grumble(ruler)`. "Never bring that fucking cretin in here again. He didn't drop the bomb. I did.
  - Representative operation: `if ~isa(ruler,'matlab.graphics.axis.decorator.NumericRuler')`.
  - Representative operation: `error('X axis must be numeric.')`.

## Parameters / inputs

- none

## Outputs

- updates the current axis system and installs an auto-updater

## Implementation structure

- Switches X axis tick labels to engineering notation by setting
- the numeric-ruler exponent to a multiple of three and installs
- a pan/zoom auto-updater. Syntax:
- kxtickfix()
- none
- updates the current axis system and installs an auto-updater
- Current X axis ruler
- Check consistency
- Preserve any existing limit-change callback
- Install the automatic exponent updater
- Set the current exponent
- Local update callback

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isappdata()`, `getappdata()`, `isfield()`, `isequal()`, `local_callback()`, `setappdata()`, `local_update()`, `ishandle()`, `iscell()`, `feval()`, `ischar()`, `isstring()`, `evalin()`, `char()`, `strcmp()`.
