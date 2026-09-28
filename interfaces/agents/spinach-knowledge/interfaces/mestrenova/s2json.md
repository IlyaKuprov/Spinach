# interfaces/mestrenova/s2json.m

- Signature: `s2json(file_name,sys,inter,parameters,fid)`

## Purpose

Writes Spinach system, interaction, parameter, and free-induction-decay data to a JSON file that can be imported by MestreNova. The function validates the inputs and calls `savejson('spinach',spinach,file_name)`.

## Inputs

- `file_name`: output file name.
- `sys`, `inter`, `parameters`: Spinach input structures.
- `fid`: a structure or matrix representing the free induction decay. For Fourier-transform-only data, provide a complex matrix. For 2D States quadrature data (for example, NOESY), provide a structure with `fid.cos` and `fid.sin` matrices.

## Output

The function writes a file; it has no returned output. The JSON root is `spinach`, with fields `sys`, `inter`, `parameters`, and `fid`.

## Source

[s2json.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=s2json.m)
