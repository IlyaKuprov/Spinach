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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Get the current axis object; implemented by `ax=gca()`.
- Lines 22-23: Remove any previous Spinach box overlay; implemented by `if isappdata(ax,'SpinachKBox')`.
- Lines 36-37: Delete orphaned overlay axes from the rejected implementation; implemented by `delete(findall(get(ax,'Parent'),'Type','axes','Tag','SpinachKBox'))`.
- Lines 39-40: Create the box line in the plot axes; implemented by `set(ax,'Box','off','Layer','top')`.
- Lines 46-47: Exclude the box line from autoscaling; implemented by `lim_props={'XLimInclude','YLimInclude','ZLimInclude'}`.
- Lines 54-55: Update only on data-space changes, not on camera motion; implemented by `props={'XLim','YLim','ZLim','LineWidth','XColor','Children'}`.
- Lines 67-68: Draw the box; implemented by `local_update(ax)`.

### Control flow inferred from the code

- Line 23: conditional branch on `isappdata(ax,'SpinachKBox')`.
- Line 25: conditional branch on `isfield(state,'listeners')`.
- Line 28: conditional branch on `isfield(state,'box_line')&&ishandle(state.box_line)`.
- Line 31: conditional branch on `isfield(state,'box_axes')&&ishandle(state.box_axes)`.
- Line 48: `for` loop over `n=1:numel(lim_props)`.
- Line 49: conditional branch on `isprop(state.box_line,lim_props{n})`.
- Line 57: `for` loop over `n=1:numel(props)`.
- Line 59: conditional branch on `~isempty(prop)&&prop.SetObservable`.

### Key state/data transformations

- Lines 20: computes `ax` using `ax=gca()`.
- Lines 24: computes `state` using `state=getappdata(ax,'SpinachKBox'); rmappdata(ax,'SpinachKBox')`.
- Lines 41-44: computes `state.box_line` using `state.box_line=line(ax,nan,nan,'Color',get(ax,'XColor'), 'LineWidth',get(ax,'LineWidth'),'Clipping','off', 'HandleVisibility','off','HitTest','off', 'PickableParts','non…`.
- Lines 47: computes `lim_props` using `lim_props={'XLimInclude','YLimInclude','ZLimInclude'}`.
- Lines 55: computes `props` using `props={'XLim','YLim','ZLim','LineWidth','XColor','Children'}`.
- Lines 56: computes `state.listeners` using `state.listeners={}`.
- Lines 58: computes `prop` using `prop=findprop(ax,props{n})`.
- Lines 60-61: computes `state.listeners{end+1}` using `state.listeners{end+1}=addlistener(ax,props{n},'PostSet', @(~,~)local_update(ax))`.
- Lines 65: computes `state.busy` using `state.busy=false; setappdata(ax,'SpinachKBox',state)`.

### Local helper functions

- Line 73: `local_update()` — `function local_update(ax)`. Ignore inactive or recursive updates
  - Representative operation: `if ~ishandle(ax)||~isappdata(ax,'SpinachKBox'), return, end`.
  - Representative operation: `state=getappdata(ax,'SpinachKBox')`.
- Line 117: `local_seg()` — `function [x,y,z]=local_seg(x,y,z,a,b)`. Add one segment
  - Representative operation: `x=[x a(1) b(1) nan]`.
  - Representative operation: `y=[y a(2) b(2) nan]`.

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
