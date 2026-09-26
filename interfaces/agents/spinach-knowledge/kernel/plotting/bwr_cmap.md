# kernel/plotting/bwr_cmap.m

- Signature: `cmap=bwr_cmap()`

## Purpose

Blue -> White -> Red colour map with 255 points and white colour corresponding to zero. Syntax: cmap=bwr_cmap() The output is 255x3 RGB column that starts at blue, goes into white and then into red in with a quadra- tic bend.

## Physical / mathematical content

## Numerical / algorithmic content

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
