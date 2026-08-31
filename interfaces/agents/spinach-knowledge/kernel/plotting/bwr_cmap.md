# kernel/plotting/bwr_cmap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/bwr_cmap.m`
- Signature: `cmap=bwr_cmap()`
- Total lines: 42

## Purpose

Blue -> White -> Red colour map with 255 points and white colour corresponding to zero. Syntax: cmap=bwr_cmap() The output is 255x3 RGB column that starts at blue, goes into white and then into red in with a quadra- tic bend.

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Preallocate the map; implemented by `cmap=zeros(255,3)`.
- Lines 23-24: Rise from blue to white; implemented by `cmap(1:128,1)=linspace(0,1,128)`.
- Lines 28-29: Rise from white to red; implemented by `cmap(128:255,1)=1`.
- Lines 33-34: Improve contrast; implemented by `cmap=cmap.^2`.

### Key state/data transformations

- Lines 21: computes `cmap` using `cmap=zeros(255,3)`.
- Lines 24: computes `cmap(1:128,1)` using `cmap(1:128,1)=linspace(0,1,128)`.
- Lines 25: computes `cmap(1:128,2)` using `cmap(1:128,2)=linspace(0,1,128)`.
- Lines 26: computes `cmap(1:128,3)` using `cmap(1:128,3)=1`.
- Lines 29: computes `cmap(128:255,1)` using `cmap(128:255,1)=1`.
- Lines 30: computes `cmap(128:255,2)` using `cmap(128:255,2)=fliplr(linspace(0,1,128))`.
- Lines 31: computes `cmap(128:255,3)` using `cmap(128:255,3)=fliplr(linspace(0,1,128))`.

## Outputs

- cmap -colour map in Matlab format

## Implementation structure

- Blue -> White -> Red colour map with 255 points and
- white colour corresponding to zero. Syntax:
- cmap=bwr_cmap()
- The output is 255x3 RGB column that starts at blue,
- goes into white and then into red in with a quadra-
- tic bend.
- cmap -colour map in Matlab format
- Preallocate the map
- Rise from blue to white
- Rise from white to red
- Improve contrast
- The worst thing I can be is the same as everybody

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cmap()`, `fliplr()`.
