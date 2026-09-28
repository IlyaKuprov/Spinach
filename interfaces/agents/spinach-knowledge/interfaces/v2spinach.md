# interfaces/v2spinach.m

- Signature: `vdata=v2spinach(inpath)`

## Purpose

Imports time-domain NMR data recorded by Varian and Agilent instruments. It reads the binary `fid` file and the `procpar` parameter file from the experiment directory.

## Parameters / inputs

- `inpath` — character string giving the path to the experiment directory containing `fid` and `procpar`.

## Outputs

- `vdata.fid` — matrix of complex free induction decays, one trace per column. Traces are ordered within each data block first, followed by the data blocks.
- `vdata.procpar` — structure containing every parameter found in the `procpar` file. Numeric parameters are column vectors; string parameters are character strings or cell arrays.
- `vdata.dirname` — experiment directory name.
- `vdata.nblocks` — number of data blocks in the `fid` file.
- `vdata.ntraces` — number of traces in each data block.
- `vdata.np` — stored points per trace, counting real and imaginary values separately.
- `vdata.ebytes` — bytes per stored data point.
- `vdata.tbytes` — bytes per trace.
- `vdata.bbytes` — bytes per data block.
- `vdata.version_id` — VnmrJ version identifier.
- `vdata.status` — `fid` file status bitfield.
- `vdata.nbheaders` — number of block headers per data block.
- `vdata.block` — block-header array with fields `scale`, `status`, `bitstatus`, `index`, `mode`, `ctcount`, `lpval`, `rpval`, `lvl`, and `tlt`.
- `vdata.npoints` — complex points per trace.
- `vdata.arraydim` — number of elements in the parameter array.
- `vdata.sfrq` — spectrometer frequency, MHz.
- `vdata.at` — acquisition time, seconds.
- `vdata.sw_ppm` — spectral width, ppm.
- `vdata.spec_start` — lower edge of the spectrum, ppm.
- `vdata.ni` — indirect-dimension increment count; present when specified in `procpar`.
- `vdata.grad_amps` — diffusion gradient amplitudes, T/m; present when `procpar` contains a gradient array and a gradient calibration factor.
- `vdata.big_delta` — diffusion delay, seconds; present when specified in `procpar`.
- `vdata.small_delta` — diffusion-encoding gradient duration, seconds; present when specified in `procpar`.
- `vdata.gamma` — magnetogyric ratio used by the diffusion sequence, rad/(s*T); present when specified in `procpar`.
- `vdata.dosy_const` — Stejskal–Tanner time factor, `gamma^2 * delta^2 * (DELTA - delta/3)`, present when specified in `procpar`. In a pulsed-field-gradient experiment, the signal attenuation is `exp(-dosy_const*g^2*D)`, where `g` is the gradient amplitude in T/m and `D` is the diffusion coefficient in m^2/s.

## Implementation overview

The function checks that `inpath` is a character string naming an existing directory containing both required files. It opens `fid` as a big-endian binary file, reads its header, and checks the trace and block sizes for consistency. It then reads each block header, skips any additional hypercomplex block headers, selects the stored numeric type from the block status bits, and assembles the interleaved real and imaginary values into complex traces.

The `procpar` reader stores numeric values as vectors and string values as character strings or cell arrays. The function then derives the frequency and spectral fields and, when the corresponding parameters are present, the indirect-dimension count, gradient amplitudes, diffusion timing, magnetogyric ratio, and Stejskal–Tanner factor.

## Provenance

Adapted from the `varianimport()` function of the GNAT package by Dr. Mathias Nilsson, School of Chemistry, University of Manchester, Oxford Road, Manchester M13 9PL, UK (mathias.nilsson@manchester.ac.uk). Contact: ilya.kuprov@weizmann.ac.il.

<https://spindynamics.org/wiki/index.php?title=v2spinach.m>
