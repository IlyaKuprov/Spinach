# interfaces/gaussian/gslice.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/gaussian/gslice.m`
- Signature: `gslice()`
- Total lines: 103

## Purpose

Slices a Gaussian geometry scan log into property calculation inputs at the energy minimum geometries. The function asks for the log and for a text file containing the property calcula- tion header. Some paths are hard-coded; edit as appropriate.

## Physical / mathematical content

- Gaussian interfaces. These parse quantum-chemistry output into spin Hamiltonian ingredients such as hyperfine, shielding, or exchange parameters.

## Numerical / algorithmic content

- The file also defines local helper function(s): `dump()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Slices a Gaussian geometry scan log into property calculation
- inputs at the energy minimum geometries. The function asks for
- the log and for a text file containing the property calcula-
- tion header. Some paths are hard-coded; edit as appropriate.
- Assign atomic symbols
- Read the log file
- Locate and the stationary point reports
- Locate standard orientation entry points
- Locate standard orientation end points
- Read the standard orientations
- Get the header
- Write the inputs

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `uigetfile()`, `textread()`, `strcmp()`, `deblank()`, `g03_output()`, `num2str()`, `std_geom_end()`, `std_geom_start()`, `cell2mat()`, `textscan()`, `fopen()`, `dump()`, `fclose()`.
