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
