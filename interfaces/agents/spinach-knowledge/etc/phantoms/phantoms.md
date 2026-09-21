# etc/phantoms/phantoms.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/phantoms/phantoms.m`
- Signature: `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)`
- Total lines: 169

## Purpose

MRI phantom library. Syntax: [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- ph_name -character string giving the name of the
- phantom (see the function text)

## Outputs

- R1Ph -a cube of R1 values
- R2Ph -a cube of R2 values
- PDPh -a cube of PD values
- dims -row vector of three cube dimensions, m
- npts -row vector of three cube dimensions, points

## Implementation structure

- MRI phantom library. Syntax:
- [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)
- ph_name -character string giving the name of the
- phantom (see the function text)
- R1Ph -a cube of R1 values
- R2Ph -a cube of R2 values
- PDPh -a cube of PD values
- dims -row vector of three cube dimensions, m
- npts -row vector of three cube dimensions, points
- Check consistency
- Get own location
- Load the breast phantom

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mfilename()`, `load()`, `PDPh()`, `R1Ph()`, `R2Ph()`, `isinf()`, `ischar()`.
