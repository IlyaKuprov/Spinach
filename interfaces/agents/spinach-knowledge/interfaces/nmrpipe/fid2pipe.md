# interfaces/nmrpipe/fid2pipe.m

- Signature: `fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)`

## Purpose

Exports phase-sensitive 2D Spinach free induction decays into native NMRPipe time-domain files. Syntax: fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)

## Physical / mathematical content

- NMRPipe interfaces. These files write or translate Spinach time-domain and spectral data into formats used by NMRPipe-style processing workflows.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system structure
- file_root -output file name root without shell metacharacters;
- the function writes <file_root>_pos.fid and
- <file_root>_neg.fid
- fid -structure with fid.pos and fid.neg complex matrices;
- rows are direct-dimension points, and columns are
- indirect-dimension points
- parameters -pulse sequence parameters structure with fields
- sweep, npoints, offset, and spins
- nmrpipe_root -NMRPipe installation root containing com/txt2pipe.tcl
- and nmrbin.linux212_64

## Outputs

- this function writes two NMRPipe files
- Note: the exported files preserve the two Spinach echo/anti-echo
- components; recombination is done in the accompanying NMRPipe
- processing script after the direct-dimension Fourier transform

## Implementation structure

- Exports phase-sensitive 2D Spinach free induction decays into native
- NMRPipe time-domain files. Syntax:
- fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)
- spin_system -Spinach spin system structure
- file_root -output file name root without shell metacharacters;
- the function writes <file_root>_pos.fid and
- <file_root>_neg.fid
- fid -structure with fid.pos and fid.neg complex matrices;
- rows are direct-dimension points, and columns are
- indirect-dimension points
- parameters -pulse sequence parameters structure with fields
- sweep, npoints, offset, and spins
