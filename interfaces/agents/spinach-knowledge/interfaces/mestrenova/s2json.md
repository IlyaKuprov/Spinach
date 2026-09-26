# interfaces/mestrenova/s2json.m

- Signature: `s2json(file_name,sys,inter,parameters,fid)`

## Purpose

Writes the parameters structure and the free induction decay into a JSON file that can be imported by MestreNova. Syntax: s2json(file_name,parameters,fid_matrices)

## Physical / mathematical content

- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
