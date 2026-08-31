# kernel/pulses/bruker_write.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/bruker_write.m`
- Signature: `bruker_write(X,Y,dt,file_name)`
- Total lines: 108

## Purpose

Saves pulses in Bruker format. The result is a text file with a list of amplitudes and phases, usable in TopSpin. Syntax: bruker_write(X,Y,dt,file_name)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(X,Y,dt,file_name)`.
- Lines 34-35: Get amplitudes and phases; implemented by `[amp,phi]=cartesian2polar(X,Y)`.
- Lines 37-38: Wrap phases and convert to degrees; implemented by `phi=rad2deg(wrapTo2Pi(phi))`.
- Lines 40-41: Interval count, min and max; implemented by `n_ints=numel(X); T=n_ints*dt`.
- Lines 45-46: Date and time; implemented by `date=datetime('now','Format','dd/MM/yyyy')`.
- Lines 49-73: Write shaped pulse lines; implemented by `lines = ["##TITLE= Optimal control pulse", "##JCAMP-DX= 5.00 Bruker JCAMP library", "##DATA TYPE= Shape Data", "##ORIGIN= SPINACH", "##OWNER= <demo>", "##DATE=" + string…`.
- Lines 75-76: Append to output file; implemented by `writelines(lines,file_name)`.

### Key state/data transformations

- Lines 35: computes `[amp,phi]` using `[amp,phi]=cartesian2polar(X,Y)`.
- Lines 38: computes `phi` using `phi=rad2deg(wrapTo2Pi(phi))`.
- Lines 41: computes `n_ints` using `n_ints=numel(X); T=n_ints*dt`.
- Lines 42: computes `min_amp` using `min_amp=min(amp); max_amp=max(amp)`.
- Lines 43: computes `min_phi` using `min_phi=min(phi); max_phi=max(phi)`.
- Lines 46: computes `date` using `date=datetime('now','Format','dd/MM/yyyy')`.
- Lines 47: computes `time` using `time=datetime('now','Format','HH:mm:ss')`.
- Lines 50-73: computes `lines` using `lines = ["##TITLE= Optimal control pulse", "##JCAMP-DX= 5.00 Bruker JCAMP library", "##DATA TYPE= Shape Data", "##ORIGIN= SPINACH", "##OWNER= <demo>", "##DATE=" + string…`.
- Lines 51-73: computes `"##JCAMP-DX` using `"##JCAMP-DX= 5.00 Bruker JCAMP library", "##DATA TYPE= Shape Data", "##ORIGIN= SPINACH", "##OWNER= <demo>", "##DATE=" + string(date), "##TIME=" + string(time), "##$SHAPE…`.
- Lines 52-73: computes `"##DATA TYPE` using `"##DATA TYPE= Shape Data", "##ORIGIN= SPINACH", "##OWNER= <demo>", "##DATE=" + string(date), "##TIME=" + string(time), "##$SHAPE_PARAMETERS=", "##MINX= " + num2str(min_a…`.
- Lines 53-73: computes `"##ORIGIN` using `"##ORIGIN= SPINACH", "##OWNER= <demo>", "##DATE=" + string(date), "##TIME=" + string(time), "##$SHAPE_PARAMETERS=", "##MINX= " + num2str(min_amp), "##MAXX= " + num2str(m…`.
- Lines 54-73: computes `"##OWNER` using `"##OWNER= <demo>", "##DATE=" + string(date), "##TIME=" + string(time), "##$SHAPE_PARAMETERS=", "##MINX= " + num2str(min_amp), "##MAXX= " + num2str(max_amp), "##MINY= " +…`.
- Lines 55-73: computes `"##DATE` using `"##DATE=" + string(date), "##TIME=" + string(time), "##$SHAPE_PARAMETERS=", "##MINX= " + num2str(min_amp), "##MAXX= " + num2str(max_amp), "##MINY= " + num2str(min_phi),…`.
- Lines 56-73: computes `"##TIME` using `"##TIME=" + string(time), "##$SHAPE_PARAMETERS=", "##MINX= " + num2str(min_amp), "##MAXX= " + num2str(max_amp), "##MINY= " + num2str(min_phi), "##MAXY= " + num2str(max_p…`.
- Lines 57-73: computes `"##$SHAPE_PARAMETERS` using `"##$SHAPE_PARAMETERS=", "##MINX= " + num2str(min_amp), "##MAXX= " + num2str(max_amp), "##MINY= " + num2str(min_phi), "##MAXY= " + num2str(max_phi), "##$SHAPE_EXMODE= Uni…`.
- Lines 58-73: computes `"##MINX` using `"##MINX= " + num2str(min_amp), "##MAXX= " + num2str(max_amp), "##MINY= " + num2str(min_phi), "##MAXY= " + num2str(max_phi), "##$SHAPE_EXMODE= Universal", "##$SHAPE_TOTRO…`.
- Lines 59-73: computes `"##MAXX` using `"##MAXX= " + num2str(max_amp), "##MINY= " + num2str(min_phi), "##MAXY= " + num2str(max_phi), "##$SHAPE_EXMODE= Universal", "##$SHAPE_TOTROT=", "##$SHAPE_TYPE= Excitation…`.
- Lines 60-73: computes `"##MINY` using `"##MINY= " + num2str(min_phi), "##MAXY= " + num2str(max_phi), "##$SHAPE_EXMODE= Universal", "##$SHAPE_TOTROT=", "##$SHAPE_TYPE= Excitation", "##$SHAPE_USER_DEF=", "##$SH…`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(X,Y,dt,file_name)`.
  - Representative operation: `if (~isnumeric(X))||(~isreal(X))||(~iscolumn(X))`.
  - Representative operation: `error('X must be a real column vector.')`.

## Parameters / inputs

- X -in-phase channel pulse amplitudes, a
- column vector in Hz
- Y -in-phase channel pulse amplitudes, a
- column vector in Hz
- dt -time slice duration, seconds
- file_name -name of output file with .txt exten-
- sion, a character string

## Outputs

- the function writes an ASCII text file

## Implementation structure

- Saves pulses in Bruker format. The result is a text file with a
- list of amplitudes and phases, usable in TopSpin. Syntax:
- bruker_write(X,Y,dt,file_name)
- X -in-phase channel pulse amplitudes, a
- column vector in Hz
- Y -in-phase channel pulse amplitudes, a
- dt -time slice duration, seconds
- file_name -name of output file with .txt exten-
- sion, a character string
- the function writes an ASCII text file
- Check consistency
- Get amplitudes and phases

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cartesian2polar()`, `rad2deg()`, `wrapTo2Pi()`, `datetime()`, `string()`, `num2str()`, `writelines()`, `writematrix()`, `iscolumn()`, `isscalar()`, `ischar()`.
