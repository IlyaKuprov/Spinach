# kernel/plotting/scale_figure.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/scale_figure.m`
- Signature: `scale_figure(by)`
- Total lines: 55

## Purpose

Scales the current figure from the default size by the factors provided by the user. Syntax: scale_figure(by)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Check consistency; implemented by `grumble(by)`.
- Lines 20-21: Get figure location; implemented by `loc=get(gcf,'Position')`.
- Lines 23-24: Get figure centroid; implemented by `old_fig_cent=[loc(1)+loc(3)/2 loc(2)+loc(4)/2]`.
- Lines 26-27: Get default figure size; implemented by `loc=get(0,'defaultfigureposition')`.
- Lines 30-31: Scale the figure; implemented by `new_fig_cent=old_fig_cent`.
- Lines 34-35: Update figure parameters; implemented by `set(gcf,'Position',[new_fig_cent-new_fig_size/2 new_fig_size])`.

### Key state/data transformations

- Lines 21: computes `loc` using `loc=get(gcf,'Position')`.
- Lines 24: computes `old_fig_cent` using `old_fig_cent=[loc(1)+loc(3)/2 loc(2)+loc(4)/2]`.
- Lines 28: computes `old_fig_size` using `old_fig_size=[loc(3) loc(4)]`.
- Lines 31: computes `new_fig_cent` using `new_fig_cent=old_fig_cent`.
- Lines 32: computes `new_fig_size` using `new_fig_size=by.*old_fig_size`.

### Local helper functions

- Line 40: `grumble()` — `function grumble(by)`. During his research visit to Southampton University to study the transport of fluorinated carbohydrates across erythrocyte membra-
  - Representative operation: `if (~isnumeric(by))||(~isreal(by))||(~isrow(by))||(numel(by)~=2)||any(by<=0)`.
  - Representative operation: `error('the argument must be a row vector with two positive real elements.')`.

## Parameters / inputs

- by -two-element row vector of scaling factors,
- format: [width height]

## Implementation structure

- Scales the current figure from the default size by the factors
- provided by the user. Syntax:
- scale_figure(by)
- by - two-element row vector of scaling factors,
- format: [width height]
- Check consistency
- Get figure location
- Get figure centroid
- Get default figure size
- Scale the figure
- Update figure parameters
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `get()`, `loc()`, `set()`, `isrow()`, `any()`.
