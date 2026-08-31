# kernel/plotting/plot_uf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/plot_uf.m`
- Signature: `plot_uf(spin_system,spectrum_uf,parameters)`
- Total lines: 184

## Purpose

Plotting utility for ultrafast constant-time 2D pulse sequences. Syntax: plot_uf(spin_system,spectrum_uf,parameters)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(spectrum_uf,parameters)`.
- Lines 45-46: Get Ta (duration of one single acquisition gradient); implemented by `Ta=parameters.deltat*parameters.npoints`.
- Lines 48-49: Get sweep width conventional dimension; implemented by `sweep_conv=1/(2*Ta)`.
- Lines 51-52: Get resolution conventional dimension; implemented by `res_conv=1/(2*Ta*parameters.nloops)`.
- Lines 54-55: Get axis in the conventional dimension; implemented by `axis_f2=(-sweep_conv/2:res_conv:sweep_conv/2-1)+parameters.offset(2)`.
- Lines 57-58: Get the magnetogyric ratio; implemented by `gamma=spin(parameters.spins{1})`.
- Lines 60-61: Get kmax (maximal k-value); implemented by `k_max=gamma*parameters.Ga*Ta/(2*pi)`.
- Lines 63-65: Get constant c (in according to Prog.Nucl.Magn.Reson.Spectrosc. 57,2010,241); implemented by `t_max=2*parameters.Te`.
- Lines 68-69: Get resolution uf dimension; implemented by `res_uf=abs(1/(constant_c*parameters.dims))`.
- Lines 71-72: Get sweep width uf dimension; implemented by `sweep_uf=abs(k_max/constant_c)`.
- Lines 74-75: Calculate the number of points of the Ga; implemented by `ga_points_number=parameters.dims*k_max`.
- Lines 78-79: Get size of the uf dimension; implemented by `uf_dim_size=size(spectrum_uf,1)`.
- Lines 81-82: A safety checking; implemented by `if (uf_dim_size) ~= (ga_points_number)`.
- Lines 86-87: Get axis in the uf dimension; implemented by `axis_f1=-sweep_uf/2+res_uf*(0:(uf_dim_size-1))+parameters.offset(1)`.
- Lines 89-90: Get the units; implemented by `switch parameters.axis_units`.
- Lines 97-98: Add offset between uf and conventional F1 axes; implemented by `axis_f1=axis_f1+parameters.offset_uf_cov`.
- Lines 109-110: Complain and bomb out; implemented by `error('unknown axis units.')`.
- Lines 114-115: Plot the spectrum; implemented by `contour(axis_f2,axis_f1,flipud(spectrum_uf))`.

### Control flow inferred from the code

- Line 82: conditional branch on `(uf_dim_size) ~= (ga_points_number)`.
- Line 90: dispatches on `parameters.axis_units`; cases `'ppm'`, `'Hz'`.

### Key state/data transformations

- Lines 46: computes `Ta` using `Ta=parameters.deltat*parameters.npoints`.
- Lines 49: computes `sweep_conv` using `sweep_conv=1/(2*Ta)`.
- Lines 52: computes `res_conv` using `res_conv=1/(2*Ta*parameters.nloops)`.
- Lines 55: computes `axis_f2` using `axis_f2=(-sweep_conv/2:res_conv:sweep_conv/2-1)+parameters.offset(2)`.
- Lines 58: computes `gamma` using `gamma=spin(parameters.spins{1})`.
- Lines 61: computes `k_max` using `k_max=gamma*parameters.Ga*Ta/(2*pi)`.
- Lines 65: computes `t_max` using `t_max=2*parameters.Te`.
- Lines 66: computes `constant_c` using `constant_c=(-2*t_max)/parameters.dims`.
- Lines 69: computes `res_uf` using `res_uf=abs(1/(constant_c*parameters.dims))`.
- Lines 72: computes `sweep_uf` using `sweep_uf=abs(k_max/constant_c)`.
- Lines 75: computes `ga_points_number` using `ga_points_number=parameters.dims*k_max`.
- Lines 79: computes `uf_dim_size` using `uf_dim_size=size(spectrum_uf,1)`.
- Lines 87: computes `axis_f1` using `axis_f1=-sweep_uf/2+res_uf*(0:(uf_dim_size-1))+parameters.offset(1)`.

### Local helper functions

- Line 127: `grumble()` — `function grumble(spectrum_uf,parameters)`.
  - Representative operation: `if (~isnumeric(spectrum_uf))||(~isreal(spectrum_uf))||(~ismatrix(spectrum_uf))`.
  - Representative operation: `error('the uf spectrum must be a real matrix.')`.

## Parameters / inputs

- spectrum_uf -a real matrix containing the 2D UF NMR
- spectrum
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.dims -sample dimension, m
- parameters.deltat -time step for the acquisition gradient, s
- parameters.npoints -number of points in the acquisition
- gradient
- parameters.Ga -amplitude of the acquisition gradient, T/m
- parameters.offset -two transmitter offsets for the
- conventional and uf dimension, Hz
- parameters.axis_units -axis units ('ppm' or 'Hz')
- parameters.offset_uf_cov -offset between chemical shifts of a MQ
- along the F1 dimension of a conventional
- and an UF spectra ('ppm' or 'Hz').
- Output:
- a figure with correct axis ticks axes in the UF and conventional
- dimension

## Implementation structure

- Plotting utility for ultrafast constant-time 2D pulse sequences. Syntax:
- plot_uf(spin_system,spectrum_uf,parameters)
- spectrum_uf -a real matrix containing the 2D UF NMR
- spectrum
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.dims -sample dimension, m
- parameters.deltat -time step for the acquisition gradient, s
- parameters.npoints -number of points in the acquisition
- gradient
- parameters.Ga -amplitude of the acquisition gradient, T/m
- parameters.offset -two transmitter offsets for the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `flipud()`, `kxlabel()`, `kylabel()`, `set()`, `ismatrix()`, `isfield()`, `isscalar()`, `ischar()`, `ismember()`.
