# kernel/pulses/read_wave.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/read_wave.m`
- Signature: `[A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)`
- Total lines: 95

## Purpose

Reads JCAMP-DX pulse waveform files (a few examples are distri- buted with Spinach, see /kernel/pulses/pk_files). Syntax: [A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(filename,npoints)`.
- Lines 39-40: Read the waveform; implemented by `P=mfilename('fullpath'); P=P(1:(end-9))`.
- Lines 45-46: Read the scaling factor; implemented by `line_by_line=textscan(wavefile,'%s','delimiter','\n')`.
- Lines 58-59: Get amplitude and phase; implemented by `A=waveform(:,1)/100; phi=pi*waveform(:,2)/180`.
- Lines 61-62: Resample amplitude and phase; implemented by `A=interp1(linspace(0,1,numel(A)),A,linspace(0,1,npoints),'pchip')`.
- Lines 65-66: Convert into Cartesian coordinates; implemented by `if nargout>3,[Cx,Cy]=polar2cartesian(A,phi); end`.

### Control flow inferred from the code

- Line 48: `for` loop over `n=1:numel(line_by_line)`.
- Line 49: conditional branch on `(numel(line_by_line{n})>17)&&`.
- Line 54: conditional branch on `(~exist('scaling_factor','var'))||isnan(scaling_factor)`.
- Line 66: conditional branch on `nargout>3,[Cx,Cy]=polar2cartesian(A,phi); end`.

### Key state/data transformations

- Lines 40: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-9))`.
- Lines 41: computes `wavefile` using `wavefile=fopen([P 'pk_files' filesep filename],'r')`.
- Lines 42: computes `waveform` using `waveform=textscan(wavefile,'%f, %f','CommentStyle','##')`.
- Lines 46: computes `line_by_line` using `line_by_line=textscan(wavefile,'%s','delimiter','\n')`.
- Lines 47: computes `fclose(wavefile); line_by_line` using `fclose(wavefile); line_by_line=line_by_line{1}`.
- Lines 51: computes `scaling_factor` using `scaling_factor=str2double(line_by_line{n}(19:end))`.
- Lines 55: computes `error('the scaling factor string '' ##$SHAPE_INTEGFAC` using `error('the scaling factor string '' ##$SHAPE_INTEGFAC= '' must be present in the pk_file.')`.
- Lines 59: computes `A` using `A=waveform(:,1)/100; phi=pi*waveform(:,2)/180`.
- Lines 63: computes `phi` using `phi=interp1(linspace(0,1,numel(phi)),phi,linspace(0,1,npoints),'pchip')`.

### Local helper functions

- Line 71: `grumble()` — `function grumble(filename,npoints)`. I condemn Christianity... It is, to me, the greatest of all imaginable corruptions. [...] To breed in humans a self-contradiction, an art of
  - Representative operation: `if ~ischar(filename)`.
  - Representative operation: `error('filename must be a character string.')`.

## Parameters / inputs

- filename -a string containing the name of the file
- npoints -waveform upsampling or downsampling is
- performed to this number of points

## Outputs

- A -polar amplitude at each slice
- phi -polar phase at each slice, radians
- Cx -Cartesian amplitude in X at each slice
- Cy -Cartesian amplitude in Y at each slice
- scaling factor -scaling factor for a given pulse shape
- Note: put your own pulses into /kernel/pulses/pk_files; please
- also consider sending them to us.

## Implementation structure

- Reads JCAMP-DX pulse waveform files (a few examples are distri-
- buted with Spinach, see /kernel/pulses/pk_files). Syntax:
- [A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)
- filename - a string containing the name of the file
- npoints - waveform upsampling or downsampling is
- performed to this number of points
- A - polar amplitude at each slice
- phi - polar phase at each slice, radians
- Cx - Cartesian amplitude in X at each slice
- Cy - Cartesian amplitude in Y at each slice
- scaling factor - scaling factor for a given pulse shape
- Note: put your own pulses into /kernel/pulses/pk_files; please

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mfilename()`, `fopen()`, `textscan()`, `cell2mat()`, `frewind()`, `fclose()`, `strcmp()`, `str2double()`, `exist()`, `isnan()`, `waveform()`, `polar2cartesian()`, `ischar()`.
