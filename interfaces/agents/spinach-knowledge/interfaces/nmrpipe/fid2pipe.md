# interfaces/nmrpipe/fid2pipe.m

- Signature: `fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)`

## Purpose

Exports phase-sensitive 2D Spinach free induction decays as native NMRPipe time-domain files. The two echo/anti-echo components are preserved; recombination is performed by the accompanying NMRPipe processing script after the direct-dimension Fourier transform.

## Inputs

- `spin_system`: Spinach spin-system structure.
- `file_root`: output name root without shell metacharacters; files are `<file_root>_pos.fid` and `<file_root>_neg.fid`.
- `fid`: structure containing complex matrices `fid.pos` and `fid.neg`; rows are direct-dimension points and columns are indirect-dimension points.
- `parameters`: pulse-sequence structure with `sweep`, `npoints`, `offset`, and `spins` fields.
- `nmrpipe_root`: NMRPipe installation root containing `com/txt2pipe.tcl` and `nmrbin.linux212_64`.

## Encoding

The indirect dimension follows the NMRPipe hypercomplex convention: two records are written per increment, so `yN` is twice `yT`.

## Output and source

The function writes two NMRPipe files and returns no value. [fid2pipe.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=fid2pipe.m)
