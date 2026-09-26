# interfaces/v2spinach.m

- Signature: `vdata=v2spinach(inpath)`

## Purpose

Imports time-domain NMR data recorded by Varian and Agilent inst- ruments: reads the binary fid file and the procpar parameter file from the experiment directory. Syntax: vdata=v2spinach(inpath)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- inpath -character string with the path to the experiment
- directory containing fid and procpar files

## Outputs

- vdata.fid -matrix of complex free induction decays,
- one per column, in the order they appear
- in the file: traces within a data block
- run first, data blocks run second
- vdata.procpar -structure with every parameter found in
- the procpar file: numeric parameters as
- column vectors, string parameters as
- character strings or cell arrays thereof
- vdata.dirname -experiment directory name
- vdata.nblocks -number of data blocks in the fid file
- vdata.ntraces -number of traces in each data block
- vdata.np -stored points per trace, real and imagi-
- nary parts counted separately
- vdata.ebytes -bytes per stored data point
- vdata.tbytes -bytes per trace
- vdata.bbytes -bytes per data block
- vdata.version_id -VnmrJ version identifier
- vdata.status -fid file status bitfield
- vdata.nbheaders -number of block headers per data block
- vdata.block -block header array with fields: scale,
- status, bitstatus, index, mode, ctcount,
- lpval, rpval, lvl, and tlt
- vdata.npoints -complex points per trace
- vdata.arraydim -number of elements in the parameter array
- vdata.sfrq -spectrometer frequency, MHz
- vdata.at -acquisition time, seconds
- vdata.sw_ppm -spectral width, ppm
- vdata.spec_start -lower edge of the spectrum, ppm
- vdata.ni -increment count of the indirect dimensi-
- on, present when procpar specifies it
- vdata.grad_amps -diffusion gradient amplitudes, T/m, pre-
- sent when procpar holds a gradient array
- and a gradient calibration factor
- vdata.big_delta -diffusion delay, seconds, present when
- procpar specifies it
- vdata.small_delta -diffusion encoding gradient duration,
- seconds, present when procpar specifies it
- vdata.gamma -magnetogyric ratio used by the diffusion
- sequence, rad/(s*T), present when procpar
- specifies it
- vdata.dosy_const -Stejskal-Tanner time factor: gamma^2 mul-
- tiplied by the effective delta^2*(DELTA-
- delta/3) term of the pulse sequence, pre-
- sent when procpar specifies it; the sig-
- nal attenuation in a pulsed field gradi-
- ent experiment is exp(-dosy_const*g^2*D)
- where g is the gradient amplitude in T/m
- and D is the diffusion coefficient
- in m^2/s
- Adapted from the varianimport() function of the GNAT package by:
- Dr. Mathias Nilsson
- School of Chemistry, University of Manchester,
- Oxford Road, Manchester M13 9PL, UK

## Implementation structure

- Imports time-domain NMR data recorded by Varian and Agilent inst-
- ruments: reads the binary fid file and the procpar parameter file
- from the experiment directory. Syntax:
- vdata=v2spinach(inpath)
- inpath - character string with the path to the experiment
- directory containing fid and procpar files
- vdata.fid -matrix of complex free induction decays,
- one per column, in the order they appear
- in the file: traces within a data block
- run first, data blocks run second
- vdata.procpar -structure with every parameter found in
- the procpar file: numeric parameters as
