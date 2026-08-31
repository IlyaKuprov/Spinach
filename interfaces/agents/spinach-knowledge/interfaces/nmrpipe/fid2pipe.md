# interfaces/nmrpipe/fid2pipe.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/nmrpipe/fid2pipe.m`
- Signature: `fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)`
- Total lines: 195

## Purpose

Exports phase-sensitive 2D Spinach free induction decays into native NMRPipe time-domain files. Syntax: fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)

## Physical / mathematical content

- NMRPipe interfaces. These files write or translate Spinach time-domain and spectral data into formats used by NMRPipe-style processing workflows.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system,file_root,fid,parameters,nmrpipe_root)`.
- Lines 41-42: Locate the NMRPipe program directories; implemented by `nmrbin=fullfile(nmrpipe_root,'nmrbin.linux212_64')`.
- Lines 46-47: Set the NMRPipe runtime environment; implemented by `setenv('PATH',[nmrbin pathsep nmrcom pathsep getenv('PATH')])`.
- Lines 57-58: Set the OpenWindows path if it is present; implemented by `if isfolder(fullfile(nmrbin,'openwin'))`.
- Lines 62-63: Convert spectrometer frequencies into MHz; implemented by `obs_f2=abs(spin(parameters.spins{2})*spin_system.inter.magnet/(2*pi))*1e-6`.
- Lines 66-67: Convert carrier offsets into signed ppm; implemented by `car_f2=hz2ppm(parameters.offset(2),spin_system.inter.magnet,parameters.spins{2})`.
- Lines 70-71: Get acquisition dimensions; implemented by `npts_f2=size(fid.pos,1)`.
- Lines 74-75: Export echo and anti-echo components separately; implemented by `branches={'pos','neg'}`.
- Lines 78-79: Select the current component; implemented by `data=fid.(branches{n})`.
- Lines 81-82: Write the text input for txt2pipe; implemented by `txt_file=[file_root '_' branches{n} '.txt']`.
- Lines 94-95: Convert the text file into NMRPipe format; implemented by `pipe_file=[file_root '_' branches{n} '.fid']`.
- Lines 112-113: Remove the temporary text file; implemented by `delete(txt_file)`.
- Lines 115-116: Report NMRPipe conversion failures; implemented by `if status~=0`.

### Control flow inferred from the code

- Line 58: conditional branch on `isfolder(fullfile(nmrbin,'openwin'))`.
- Line 76: `for` loop over `n=1:numel(branches)`.
- Line 84: conditional branch on `file_id<0`.
- Line 87: `for` loop over `k=1:npts_f1`.
- Line 88: `for` loop over `m=1:npts_f2`.
- Line 116: conditional branch on `status~=0`.

### Key state/data transformations

- Lines 42: computes `nmrbin` using `nmrbin=fullfile(nmrpipe_root,'nmrbin.linux212_64')`.
- Lines 43: computes `nmrcom` using `nmrcom=fullfile(nmrpipe_root,'com')`.
- Lines 44: computes `nmrtcl` using `nmrtcl=fullfile(nmrpipe_root,'nmrtcl')`.
- Lines 63: computes `obs_f2` using `obs_f2=abs(spin(parameters.spins{2})*spin_system.inter.magnet/(2*pi))*1e-6`.
- Lines 64: computes `obs_f1` using `obs_f1=abs(spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi))*1e-6`.
- Lines 67: computes `car_f2` using `car_f2=hz2ppm(parameters.offset(2),spin_system.inter.magnet,parameters.spins{2})`.
- Lines 68: computes `car_f1` using `car_f1=hz2ppm(parameters.offset(1),spin_system.inter.magnet,parameters.spins{1})`.
- Lines 71: computes `npts_f2` using `npts_f2=size(fid.pos,1)`.
- Lines 72: computes `npts_f1` using `npts_f1=size(fid.pos,2)`.
- Lines 75: computes `branches` using `branches={'pos','neg'}`.
- Lines 79: computes `data` using `data=fid.(branches{n})`.
- Lines 82: computes `txt_file` using `txt_file=[file_root '_' branches{n} '.txt']`.
- Lines 83: computes `file_id` using `file_id=fopen(txt_file,'w')`.
- Lines 95: computes `pipe_file` using `pipe_file=[file_root '_' branches{n} '.fid']`.
- Lines 96-109: computes `command` using `command=[fullfile(nmrcom,'txt2pipe.tcl') ' -in ' txt_file ' -xy -complex -time -xN ' num2str(2*npts_f2) ' -xT ' num2str(npts_f2) ' -xMODE Complex' ' -xSW ' num2str(param…`.
- Lines 110: computes `[status,cmdout]` using `[status,cmdout]=system(command)`.

### Local helper functions

- Line 125: `grumble()` — `function grumble(spin_system,file_root,fid,parameters,nmrpipe_root)`.
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a Spinach spin system structure.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fullfile()`, `setenv()`, `getenv()`, `isfolder()`, `spin()`, `hz2ppm()`, `fopen()`, `data()`, `fclose()`, `num2str()`, `system()`, `delete()`, `isstruct()`, `isfield()`, `iscell()`.
