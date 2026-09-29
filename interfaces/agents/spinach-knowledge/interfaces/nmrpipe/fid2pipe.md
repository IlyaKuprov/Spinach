# interfaces/nmrpipe/fid2pipe.m

- Signature: `fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)` (all five inputs are required; there are no optional arguments).

## Purpose and data contract

Exports the two phase-sensitive 2-D components in `fid.pos` and `fid.neg` as `<file_root>_pos.fid` and `<file_root>_neg.fid` NMRPipe time-domain files. The function returns no MATLAB output. It forms the indirect-dimension hypercomplex records expected for States acquisition and invokes the NMRPipe `txt2pipe.tcl` converter; it does not return a processing script.

- `spin_system`: structure with a finite, nonzero, real scalar `inter.magnet` and cell-array `comp.isotopes`.
- `file_root`: nonempty character string restricted to letters, digits, `/`, `_`, `.`, and `-`; whitespace and other shell metacharacters are rejected. Any directory in the root must already exist.
- `fid`: structure with finite numeric matrix fields `pos` and `neg` of identical size. The documented use is complex echo/anti-echo FIDs; the runtime guard checks numeric finite matrices and matching shape but does not require complex values.
- `parameters`: structure with `sweep`, `npoints`, `offset`, and `spins`. The first entry of these paired acquisition fields is indirect dimension F1 and the second is direct dimension F2. `sweep` must contain two positive finite real values; they are passed unchanged to the NMRPipe F1/F2 sweep-width fields. `offset` must contain two finite real values and is converted from Hz to carrier ppm with `hz2ppm`. `npoints` must contain two positive integers, with entry 1 equal to the number of columns and entry 2 equal to the number of rows in each FID matrix. `spins` must be a two-cell array of isotope names present in `spin_system.comp.isotopes`.
- `nmrpipe_root`: existing NMRPipe installation root, restricted to the same path-safe character set. It must contain `com/txt2pipe.tcl` and `nmrbin.linux212_64/nmrPipe`.

## File encoding and units

Rows of the input matrices are direct-dimension points and columns are indirect-dimension points. The data are written with the direct index, indirect record index, real part, and imaginary part. For each indirect increment, the indirect quadrature partner is formed with a `-i` multiplier for `pos` and `+i` for `neg`. NMRPipe metadata declares complex acquisition in both dimensions: direct `xN=2*npoints(2)`, `xT=npoints(2)`; indirect `yN=2*npoints(1)`, `yT=npoints(1)`, with `yMODE=Complex` and `aq2D=States`. The labels use `spins{2}` for F2 and `spins{1}` for F1. Observation frequencies are calculated from isotope gyromagnetic ratios and `inter.magnet` in MHz; carrier offsets are converted to ppm. The sweep-width values are passed through without conversion.

The implementation creates and removes a temporary text file for each component, and raises an error if opening that file or the NMRPipe conversion fails. It does not validate the completed output files after a successful converter exit.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/nmrpipe/fid2pipe.m)
- [Spin Dynamics Wiki: fid2pipe.m](https://spindynamics.org/wiki/index.php?title=fid2pipe.m)
