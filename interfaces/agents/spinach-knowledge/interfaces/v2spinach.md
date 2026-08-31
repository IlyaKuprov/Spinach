# interfaces/v2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/v2spinach.m`
- Signature: `vdata=v2spinach(inpath)`
- Total lines: 329

## Purpose

Imports time-domain NMR data recorded by Varian and Agilent inst- ruments: reads the binary fid file and the procpar parameter file from the experiment directory. Syntax: vdata=v2spinach(inpath)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 101-102: Check consistency; implemented by `grumble(inpath)`.
- Lines 104-105: Store the directory name; implemented by `vdata.dirname=inpath`.
- Lines 107-108: Open the fid file, big-endian; implemented by `fid_file=fopen([inpath filesep 'fid'],'r','b')`.
- Lines 110-111: Read the fid file header; implemented by `vdata.nblocks=fread(fid_file,1,'int32')`.
- Lines 121-123: Check the header for internal consistency; implemented by `if (vdata.tbytes~=vdata.np*vdata.ebytes)|| (vdata.bbytes~=28*vdata.nbheaders+vdata.ntraces*vdata.tbytes)`.
- Lines 127-128: Preallocate the fid array; implemented by `vdata.fid=zeros(vdata.np/2,vdata.nblocks*vdata.ntraces,'like',1i)`.
- Lines 130-133: Preallocate the block header array; implemented by `vdata.block(vdata.nblocks)=struct('scale',[],'status',[],'bitstatus',[], 'index',[],'mode',[],'ctcount',[], 'lpval',[],'rpval',[],'lvl',[],'tlt',[])`.
- Lines 135-136: Loop over the data blocks; implemented by `for n=1:vdata.nblocks`.
- Lines 138-139: Read the block header; implemented by `vdata.block(n).scale=fread(fid_file,1,'int16')`.
- Lines 150-151: Skip hypercomplex block headers; implemented by `fseek(fid_file,28*(vdata.nbheaders-1),'cof')`.
- Lines 153-154: Get the data type from the block status bits; implemented by `if vdata.block(n).bitstatus(4)==1`.
- Lines 162-163: Read the data points of the current block; implemented by `points=fread(fid_file,vdata.np*vdata.ntraces,data_type)`.
- Lines 165-166: Make sure the block was complete; implemented by `if numel(points)<vdata.np*vdata.ntraces`.
- Lines 170-171: Assemble complex fids, one trace per column; implemented by `points=reshape(points,[vdata.np vdata.ntraces])`.
- Lines 177-178: Close the fid file; implemented by `fclose(fid_file)`.
- Lines 180-181: Open the procpar file; implemented by `par_file=fopen([inpath filesep 'procpar'],'rt')`.
- Lines 183-184: Loop over the parameters; implemented by `while true`.
- Lines 186-187: Read the parameter description line; implemented by `name_line=fgetl(par_file)`.

### Control flow inferred from the code

- Line 122: conditional branch on `(vdata.tbytes~=vdata.np*vdata.ebytes)||`.
- Line 136: `for` loop over `n=1:vdata.nblocks`.
- Line 154: conditional branch on `vdata.block(n).bitstatus(4)==1`.
- Line 166: conditional branch on `numel(points)<vdata.np*vdata.ntraces`.
- Line 184: `while` loop over `true`.
- Line 190: conditional branch on `~ischar(name_line), break; end`.
- Line 193: conditional branch on `isempty(strtrim(name_line)), continue; end`.
- Line 203: conditional branch on `basic_type==1`.
- Line 210: `while` loop over `numel(par_vals)<val_count`.
- Line 225: `for` loop over `k=1:val_count`.
- Line 228: conditional branch on `k>1, val_line=fgetl(par_file); end`.
- Line 231: `while` loop over `nnz(val_line=='"')<2`.
- Line 233: conditional branch on `~ischar(next_line)`.
- Line 246: conditional branch on `val_count==1`.

### Key state/data transformations

- Lines 105: computes `vdata.dirname` using `vdata.dirname=inpath`.
- Lines 108: computes `fid_file` using `fid_file=fopen([inpath filesep 'fid'],'r','b')`.
- Lines 111: computes `vdata.nblocks` using `vdata.nblocks=fread(fid_file,1,'int32')`.
- Lines 112: computes `vdata.ntraces` using `vdata.ntraces=fread(fid_file,1,'int32')`.
- Lines 113: computes `vdata.np` using `vdata.np=fread(fid_file,1,'int32')`.
- Lines 114: computes `vdata.ebytes` using `vdata.ebytes=fread(fid_file,1,'int32')`.
- Lines 115: computes `vdata.tbytes` using `vdata.tbytes=fread(fid_file,1,'int32')`.
- Lines 116: computes `vdata.bbytes` using `vdata.bbytes=fread(fid_file,1,'int32')`.
- Lines 117: computes `vdata.version_id` using `vdata.version_id=fread(fid_file,1,'int16')`.
- Lines 118: computes `vdata.status` using `vdata.status=fread(fid_file,1,'int16')`.
- Lines 119: computes `vdata.nbheaders` using `vdata.nbheaders=fread(fid_file,1,'int32')`.
- Lines 128: computes `vdata.fid` using `vdata.fid=zeros(vdata.np/2,vdata.nblocks*vdata.ntraces,'like',1i)`.
- Lines 131-133: computes `vdata.block(vdata.nblocks)` using `vdata.block(vdata.nblocks)=struct('scale',[],'status',[],'bitstatus',[], 'index',[],'mode',[],'ctcount',[], 'lpval',[],'rpval',[],'lvl',[],'tlt',[])`.
- Lines 139: computes `vdata.block(n).scale` using `vdata.block(n).scale=fread(fid_file,1,'int16')`.
- Lines 140: computes `vdata.block(n).status` using `vdata.block(n).status=fread(fid_file,1,'int16')`.
- Lines 141: computes `vdata.block(n).bitstatus` using `vdata.block(n).bitstatus=bitget(uint16(vdata.block(n).status),1:16)`.
- Lines 142: computes `vdata.block(n).index` using `vdata.block(n).index=fread(fid_file,1,'int16')`.
- Lines 143: computes `vdata.block(n).mode` using `vdata.block(n).mode=fread(fid_file,1,'int16')`.

### Local helper functions

- Line 310: `grumble()` — `function grumble(inpath)`.
  - Representative operation: `if ~ischar(inpath)`.
  - Representative operation: `error('inpath must be a character string.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `fread()`, `bitget()`, `uint16()`, `fseek()`, `points()`, `fclose()`, `fgetl()`, `ischar()`, `strtrim()`, `textscan()`, `sscanf()`, `numbers()`, `nnz()`, `val_line()`.
