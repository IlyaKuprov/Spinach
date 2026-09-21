# kernel/utilities/expdrop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/expdrop.m`
- Signature: `drop=expdrop(from,to,duration,npoints,drop_rate)`
- Total lines: 67

## Purpose

Exponential drop function. Produces an exponential fall-off from a specified value to a specified value with the specified rate and the number of points. Syntax: drop=expdrop(from,to,duration,npoints,drop_rate)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- from -the value to drop from
- to -the value to drop to
- duration -drop duration, seconds
- npoints -the number of discretisation points
- in the drop
- drop_rate -exponential drop rate, Hz

## Outputs

- drop -a row vector with the fall-off

## Implementation structure

- Exponential drop function. Produces an exponential fall-off
- from a specified value to a specified value with the specified
- rate and the number of points. Syntax:
- drop=expdrop(from,to,duration,npoints,drop_rate)
- from -the value to drop from
- to -the value to drop to
- duration -drop duration, seconds
- npoints -the number of discretisation points
- in the drop
- drop_rate -exponential drop rate, Hz
- drop -a row vector with the fall-off
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
