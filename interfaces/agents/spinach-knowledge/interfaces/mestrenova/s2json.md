# interfaces/mestrenova/s2json.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/mestrenova/s2json.m`
- Signature: `s2json(file_name,sys,inter,parameters,fid)`
- Total lines: 70

## Purpose

Writes the parameters structure and the free induction decay into a JSON file that can be imported by MestreNova. Syntax: s2json(file_name,parameters,fid_matrices)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(file_name,sys,inter,parameters,fid)`.
- Lines 37-38: Form the output structure; implemented by `spinach.sys=sys; spinach.inter=inter`.
- Lines 42-43: Save a JSON object; implemented by `savejson('spinach',spinach,file_name)`.

### Key state/data transformations

- Lines 38: computes `spinach.sys` using `spinach.sys=sys; spinach.inter=inter`.
- Lines 39: computes `spinach.parameters` using `spinach.parameters=parameters`.
- Lines 40: computes `spinach.fid` using `spinach.fid=fid`.

### Local helper functions

- Line 48: `grumble()` — `function grumble(file_name,sys,inter,parameters,fid)`.
  - Representative operation: `if ~ischar(file_name)`.
  - Representative operation: `error('file_name must be a character string.')`.

## Parameters / inputs

- file_name -output file name
- sys -Spinach input data structure
- inter -Spinach input data structure
- parameters -Spinach input data structure
- fid -a structure or a matrix representing
- the free induction decay

## Outputs

- this function writes a file
- Notes: for data that only requires a Fourier transform, fid
- must be a complex matrix. For 2D States quadrature
- data (e.g. NOESY), fid.cos and fid.sin matrices must
- be supplied as a structure.

## Implementation structure

- Writes the parameters structure and the free induction decay
- into a JSON file that can be imported by MestreNova. Syntax:
- s2json(file_name,parameters,fid_matrices)
- file_name -output file name
- sys -Spinach input data structure
- inter -Spinach input data structure
- parameters -Spinach input data structure
- fid -a structure or a matrix representing
- the free induction decay
- this function writes a file
- must be a complex matrix. For 2D States quadrature
- data (e.g. NOESY), fid.cos and fid.sin matrices must

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `savejson()`, `ischar()`, `isstruct()`.
