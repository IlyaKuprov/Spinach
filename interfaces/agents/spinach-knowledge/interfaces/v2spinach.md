# interfaces/v2spinach.m

Source: [interfaces/v2spinach.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/v2spinach.m)

Signature: `vdata = v2spinach(inpath)`

Imports Varian/Agilent time-domain data from an experiment directory. It reads the binary fid file and text procpar file.

## Input and file handling

inpath must be a character string naming a directory that contains both fid and procpar. The implementation joins that directory with those filenames, opens fid as big-endian binary data and procpar as text, then closes both files. It is a local filesystem reader: the function contains no remote-transfer API or cache layer.

The FID header supplies the number of blocks, traces per block, stored real/imaginary points per trace, byte sizes, VnmrJ version ID, status, and block-header count. The function checks the recorded trace and block byte counts for consistency. It chooses the point encoding from each block's status flags (float32, int32, or int16), reshapes each block as stored points by traces, and interleaves real and imaginary values into complex data.

## Outputs

- vdata.fid is a complex matrix of size (np/2)×(nblocks*ntraces). Each column is one trace; trace columns from each block are appended in block order.
- vdata.dirname is the supplied directory path. Header fields are retained as vdata.nblocks, ntraces, np, ebytes, tbytes, bbytes, version_id, status, and nbheaders.
- vdata.block is a block-header struct array with scale, status, bitstatus, index, mode, ctcount, lpval, rpval, lvl, and tlt fields. vdata.npoints = np/2 gives complex points per trace.
- vdata.procpar stores every parsed parameter as a structure field: numeric parameters are numeric arrays; a single text value is a character array; multiple text values are a cell array.
- vdata.arraydim is copied from procpar.arraydim; vdata.sfrq is spectrometer frequency in MHz; vdata.at is acquisition time in seconds; vdata.sw_ppm = sw/sfrq; and vdata.spec_start = (rfp-rfl)/sfrq, in ppm.
- If present, vdata.ni is copied from procpar.ni. Gradient amplitudes are returned as vdata.grad_amps in T/m when gzlvl1 and either gcal_ or DAC_to_G are available. vdata.big_delta and vdata.small_delta copy del and gt1, respectively, in seconds. When available, vdata.gamma copies dosygamma in rad/(s*T); if dosytimecubed is also present, vdata.dosy_const = dosygamma^2 * dosytimecubed.

## Provenance and reference

Adapted from the varianimport() function of the GNAT package by Dr. Mathias Nilsson, School of Chemistry, University of Manchester, Oxford Road, Manchester M13 9PL, UK (mathias.nilsson@manchester.ac.uk). Contact: ilya.kuprov@weizmann.ac.uk.

[Spin Dynamics Wiki: v2spinach.m](https://spindynamics.org/wiki/index.php?title=v2spinach.m)
