# interfaces/gaussian/gslice.m

- Signature: `gslice()`

## Purpose

Slices a Gaussian geometry scan log into property calculation inputs at the energy minimum geometries. The function asks for the log and for a text file containing the property calcula- tion header. Some paths are hard-coded; edit as appropriate.

## Physical / mathematical content

- Gaussian interfaces. These parse quantum-chemistry output into spin Hamiltonian ingredients such as hyperfine, shielding, or exchange parameters.

## Numerical / algorithmic content

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
