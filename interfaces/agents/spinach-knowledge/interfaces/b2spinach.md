# interfaces/b2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/b2spinach.m`
- Signature: `bdata=b2spinach(inpath)`
- Total lines: 464

## Purpose

Imports time-domain NMR data recorded by Bruker instruments: reads the binary fid or ser file together with the acquisiti- on and processing parameter files from the numbered experi- ment directory. Syntax: bdata=b2spinach(inpath)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `read_jcamp()`, `read_list()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 105-106: Check consistency; implemented by `grumble(inpath)`.
- Lines 108-109: Store the directory name; implemented by `bdata.dirname=inpath`.
- Lines 111-112: Refuse data with more than four dimensions; implemented by `if isfile([inpath filesep 'acqu5s'])`.
- Lines 116-117: Detect data dimension from the status parameter files; implemented by `if isfile([inpath filesep 'acqu4s'])`.
- Lines 127-128: Read acquisition status parameters for every dimension; implemented by `bdata.acqus=read_jcamp([inpath filesep 'acqus'])`.
- Lines 139-140: Read processing parameters when present; implemented by `if isfile([inpath filesep 'pdata' filesep '1' filesep 'procs'])`.
- Lines 144-145: Count the fids in the data set; implemented by `bdata.arraydim=1`.
- Lines 150-151: Complex points per fid; implemented by `bdata.npoints=bdata.acqus.TD/2`.
- Lines 153-154: Pulse programme name and observe nucleus; implemented by `bdata.pulprog=bdata.acqus.PULPROG`.
- Lines 157-158: Magnetogyric ratio of the observe nucleus; implemented by `bdata.gamma=spin(bdata.nucleus)`.
- Lines 160-161: Spectrometer frequency, MHz; implemented by `bdata.sfrq=bdata.acqus.SFO1`.
- Lines 163-164: Spectral width, ppm; implemented by `bdata.sw_ppm=bdata.acqus.SW`.
- Lines 166-167: Acquisition time, seconds; implemented by `bdata.at=bdata.npoints/(bdata.sw_ppm*bdata.sfrq)`.
- Lines 169-170: Spectrum lower edge from processing or acquisition parameters; implemented by `if isfield(bdata,'procs')`.
- Lines 177-178: Digital filter group delay in complex points; implemented by `if isfield(bdata.acqus,'DIGMOD')&&(bdata.acqus.DIGMOD==0)`.
- Lines 180-181: Analogue filtering has no group delay; implemented by `bdata.digshift=0`.
- Lines 185-186: Modern hardware reports the group delay directly; implemented by `bdata.digshift=bdata.acqus.GRPDLY`.
- Lines 190-192: Decimation factors of older DSP firmware; implemented by `decims=[2 3 4 6 8 12 16 24 32 48 64 96 128 192 256 384 512 768 1024 1536 2048]`.

### Control flow inferred from the code

- Line 112: conditional branch on `isfile([inpath filesep 'acqu5s'])`.
- Line 117: conditional branch on `isfile([inpath filesep 'acqu4s'])`.
- Line 129: conditional branch on `bdata.ndims_data>=2`.
- Line 132: conditional branch on `bdata.ndims_data>=3`.
- Line 135: conditional branch on `bdata.ndims_data>=4`.
- Line 140: conditional branch on `isfile([inpath filesep 'pdata' filesep '1' filesep 'procs'])`.
- Line 146: conditional branch on `bdata.ndims_data>=2, bdata.arraydim=bdata.arraydim*bdata.acqu2s.TD; end`.
- Line 147: conditional branch on `bdata.ndims_data>=3, bdata.arraydim=bdata.arraydim*bdata.acqu3s.TD; end`.
- Line 148: conditional branch on `bdata.ndims_data>=4, bdata.arraydim=bdata.arraydim*bdata.acqu4s.TD; end`.
- Line 170: conditional branch on `isfield(bdata,'procs')`.
- Line 178: conditional branch on `isfield(bdata.acqus,'DIGMOD')&&(bdata.acqus.DIGMOD==0)`.
- Line 208: conditional branch on `~isfield(bdata.acqus,'DSPFVS')`.
- Line 214: conditional branch on `isempty(decim_row)`.
- Line 219: conditional branch on `(bdata.acqus.DSPFVS>=10)&&(bdata.acqus.DSPFVS<=12)`.

### Key state/data transformations

- Lines 109: computes `bdata.dirname` using `bdata.dirname=inpath`.
- Lines 118: computes `bdata.ndims_data` using `bdata.ndims_data=4`.
- Lines 128: computes `bdata.acqus` using `bdata.acqus=read_jcamp([inpath filesep 'acqus'])`.
- Lines 130: computes `bdata.acqu2s` using `bdata.acqu2s=read_jcamp([inpath filesep 'acqu2s'])`.
- Lines 133: computes `bdata.acqu3s` using `bdata.acqu3s=read_jcamp([inpath filesep 'acqu3s'])`.
- Lines 136: computes `bdata.acqu4s` using `bdata.acqu4s=read_jcamp([inpath filesep 'acqu4s'])`.
- Lines 141: computes `bdata.procs` using `bdata.procs=read_jcamp([inpath filesep 'pdata' filesep '1' filesep 'procs'])`.
- Lines 145: computes `bdata.arraydim` using `bdata.arraydim=1`.
- Lines 151: computes `bdata.npoints` using `bdata.npoints=bdata.acqus.TD/2`.
- Lines 154: computes `bdata.pulprog` using `bdata.pulprog=bdata.acqus.PULPROG`.
- Lines 155: computes `bdata.nucleus` using `bdata.nucleus=bdata.acqus.NUC1`.
- Lines 158: computes `bdata.gamma` using `bdata.gamma=spin(bdata.nucleus)`.
- Lines 161: computes `bdata.sfrq` using `bdata.sfrq=bdata.acqus.SFO1`.
- Lines 164: computes `bdata.sw_ppm` using `bdata.sw_ppm=bdata.acqus.SW`.
- Lines 167: computes `bdata.at` using `bdata.at=bdata.npoints/(bdata.sw_ppm*bdata.sfrq)`.
- Lines 171: computes `bdata.spec_start` using `bdata.spec_start=bdata.procs.OFFSET-bdata.sw_ppm`.
- Lines 181: computes `bdata.digshift` using `bdata.digshift=0`.
- Lines 191-192: computes `decims` using `decims=[2 3 4 6 8 12 16 24 32 48 64 96 128 192 256 384 512 768 1024 1536 2048]`.

### Local helper functions

- Line 336: `read_jcamp()` — `function params=read_jcamp(file_path)`. Read the file as a cell array of lines
  - Representative operation: `all_lines=regexp(fileread(file_path),'\r?\n','split')`.
  - Representative operation: `params=struct(); line_num=1`.
- Line 406: `read_list()` — `function vals=read_list(file_path)`. Read the file as a cell array of lines
  - Representative operation: `all_lines=regexp(fileread(file_path),'\r?\n','split')`.
  - Representative operation: `all_lines=strtrim(all_lines)`.
- Line 449: `grumble()` — `function grumble(inpath)`.
  - Representative operation: `if ~ischar(inpath)`.
  - Representative operation: `error('inpath must be a character string.')`.

## Parameters / inputs

- inpath -character string with the path to the numbered
- Bruker experiment directory containing the ac-
- qus file and the fid or ser file

## Outputs

- bdata.fid -matrix of complex free induction de-
- cays, one per column, in the order
- they are stored in the file; for data
- sets with three and more dimensions
- the loop order over the indirect di-
- mensions is determined by the AQSEQ
- parameter of the acqus structure
- bdata.acqus -structure with every parameter found
- in the acqus file: numeric parameters
- as scalars or column vectors, string
- parameters as character strings or
- cell arrays thereof
- bdata.acqu2s -same for the acqu2s file of data sets
- with two or more dimensions
- bdata.acqu3s -same for the acqu3s file of data sets
- with three or more dimensions
- bdata.acqu4s -same for the acqu4s file of four-di-
- mensional data sets
- bdata.procs -same for the procs file of the first
- processed data directory when present
- bdata.dirname -experiment directory name
- bdata.ndims_data -number of dimensions in the data set
- bdata.arraydim -number of fids the status parameter
- files declare as acquired
- bdata.fids_in_file -number of fid slots held by the bina-
- ry file, and the column count of the
- fid matrix; exceeds arraydim when the
- acquisition was preallocated or inter-
- rupted, and falls short of it for non-
- uniformly sampled data sets, where the
- fids follow the acquisition order of
- the sampling schedule in nus_list
- bdata.npoints -complex points per fid
- bdata.pulprog -pulse programme name
- bdata.nucleus -observe nucleus, e.g. '1H'
- bdata.gamma -magnetogyric ratio of the observe
- nucleus, rad/(s*T)
- bdata.sfrq -spectrometer frequency, MHz
- bdata.at -acquisition time, seconds
- bdata.sw_ppm -spectral width, ppm
- bdata.spec_start -lower edge of the spectrum, ppm
- bdata.digshift -group delay of the digital filter in
- complex points; the first round(dig-
- shift) points of each fid precede the
- start of the true signal
- bdata.grad_amps -gradient amplitudes, T/m, present
- when the difflist file exists
- bdata.vd_list -variable delay list, seconds, present
- when the vdlist file exists
- bdata.vc_list -variable counter list, present when
- the vclist file exists
- bdata.nus_list -non-uniform sampling schedule, pre-
- sent when the nuslist file exists
- Adapted from the brukerimport() function of the GNAT package by:
- Dr. Mathias Nilsson
- School of Chemistry, University of Manchester,
- Oxford Road, Manchester M13 9PL, UK

## Implementation structure

- Imports time-domain NMR data recorded by Bruker instruments:
- reads the binary fid or ser file together with the acquisiti-
- on and processing parameter files from the numbered experi-
- ment directory. Syntax:
- bdata=b2spinach(inpath)
- inpath - character string with the path to the numbered
- Bruker experiment directory containing the ac-
- qus file and the fid or ser file
- bdata.fid -matrix of complex free induction de-
- cays, one per column, in the order
- they are stored in the file; for data
- sets with three and more dimensions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfile()`, `read_jcamp()`, `spin()`, `isfield()`, `delays()`, `isnan()`, `fopen()`, `fread()`, `fclose()`, `data_pts()`, `read_list()`, `regexp()`, `fileread()`, `strncmp()`, `strfind()`.
