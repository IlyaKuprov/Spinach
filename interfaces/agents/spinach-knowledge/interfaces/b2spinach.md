# interfaces/b2spinach.m

- Signature: `bdata=b2spinach(inpath)`
- Return value: a structure containing the complex FID matrix, parsed acquisition metadata, and available auxiliary lists.

## Purpose and accepted input

Imports time-domain NMR data from a numbered Bruker experiment directory. `inpath` must be a character array naming an existing folder with an `acqus` file. A one-dimensional acquisition is read from `fid`; two or more dimensions use `ser`, with dimensionality inferred from the highest present `acqu2s`, `acqu3s`, or `acqu4s` status file. An `acqu5s` file causes an explicit error: more than four dimensions are unsupported.

## Returned structure

- `bdata.fid` is a complex matrix with one FID per column, in binary-file order. Each column contains `bdata.npoints` complex points. The source documents the indirect-dimension loop order for three-or-more-dimensional data as determined by `AQSEQ` in `acqus`; the returned data remain a column matrix.
- `bdata.acqus` contains parameters parsed from `acqus`; `bdata.acqu2s`, `acqu3s`, and `acqu4s` are included for the dimensions present. Numeric parameters are stored as scalars or column vectors; string values are character arrays or cells of character arrays.
- `bdata.procs` contains parameters from `pdata/1/procs` when that file exists.
- `bdata.dirname` is the input directory string; `ndims_data` is the detected dimensionality; `arraydim` is the number of FIDs declared by the indirect-dimension status files; and `fids_in_file` is the number of FID slots actually decoded from the binary file. The binary-file count can exceed the declared count when storage was preallocated or an acquisition was interrupted, and can fall short for non-uniformly sampled data.
- Acquisition metadata includes `npoints = acqus.TD/2`, pulse-program name `pulprog`, observe-nucleus label `nucleus` (for example, `'1H'`), `gamma = spin(nucleus)` in radians per second per tesla, `sfrq = acqus.SFO1` in MHz, `sw_ppm = acqus.SW` in ppm, and acquisition time `at = npoints/(sw_ppm*sfrq)` in seconds. `spec_start` is `procs.OFFSET-sw_ppm` when processing parameters are present; otherwise it is `1e6*(SFO1/BF1-1)-0.5*sw_ppm*SFO1/BF1` from the acquisition parameters.
- `digshift` records the digital-filter group delay in complex points. The importer reports this value but does not shift or trim the FIDs to apply a correction. `grad_amps` is included when `difflist` or `difflist.txt` exists and is multiplied by `0.01` to report tesla per metre. Optional `vd_list`, `vc_list`, and `nus_list` fields are read from their same-named files; the NUS schedule is documented as one row per sampled FID and one column per indirect dimension.

## Import transformations

Acquisition parameters are decoded by a local JCAMP reader. It considers private `##$` entries, removes inline `$$` comments, joins multi-line arrays marked with a `(start..end)` range, stores angle-bracket strings as character data or cell arrays, and converts fully numeric values to numeric arrays; other scalar values remain text. It reads `acqus` and the dimension-specific `acquNs` files, plus only the first processed data directory's `procs` file when present.

Binary byte order follows `BYTORDA` (`0` little-endian, `1` big-endian). `DTYPA` selects signed 32-bit integer, 32-bit float, or 64-bit float storage (`0`, `1`, or `2`, respectively). The code accounts for Bruker 1024-byte block padding, while also accepting unpadded FID strides when the file length matches. It reshapes the interleaved real values by FID, drops block padding, forms each complex point as real minus `1i` times imaginary, and applies the binary scale factor `2^NC`.

Digital-filter delay handling sets `digshift` to zero for analogue mode (`DIGMOD == 0`) and uses a reported `GRPDLY` when it is present and not `-1`. Otherwise, for non-unit `DECIM`, it uses the built-in decimation table for `DSPFVS` 10–12 or the source's closed-form expression for `DSPFVS == 13`. Missing or unsupported firmware/decimation combinations raise errors; data without digital filtering receive zero delay.

The list reader requires non-empty, whitespace-separated rows with a consistent number of entries. A single terminal `n`, `u`, `m`, or `s` is interpreted as a multiplier of `1e-9`, `1e-6`, `1e-3`, or `1`, respectively. Invalid numeric entries raise an error. The gradient list then receives the additional `0.01` scaling described above.

## Guardrails and dependencies

The importer errors for a missing directory or `acqus`, a missing `fid`/`ser`, unsupported dimensions, unsupported byte-order or data-type codes, malformed list data, or binary lengths inconsistent with the acquisition point count and block stride. It also errors if the binary file contains fewer FIDs than `arraydim`, except when a `nuslist` file is present; that exception does not itself validate schedule-to-column agreement. It relies on Spinach's `spin` function to obtain `gamma` from the observe nucleus.

## Links

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/b2spinach.m
- [Spinach Wiki: b2spinach.m](https://spindynamics.org/wiki/index.php?title=b2spinach.m)
