# kernel/optimcon/ctrl_trajan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ctrl_trajan.m`
- Signature: `ctrl_trajan(spin_system,waveform,traj_data,fidelities)`
- Total lines: 791

## Purpose

Diagnostic plotting function for optimal control module. Plots trajectory and control pulse analysis. Syntax: ctrl_trajan(spin_system,waveform,trajectory,fidelities)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: If no plotting is needed, exit; implemented by `if isempty(spin_system.control.plotting), return; end`.
- Lines 34-35: Check consistency; implemented by `grumble(spin_system,waveform,traj_data,fidelities)`.
- Lines 37-38: Trajectory is not plotted, that key exists to trigger a trajectory return; implemented by `spin_system.control.plotting=setdiff(spin_system.control.plotting,{'trajectory'})`.
- Lines 40-41: Count the plots; implemented by `n_plots=numel(spin_system.control.plotting)`.
- Lines 43-44: Exit if nothing to plot; implemented by `if n_plots==0, return; end`.
- Lines 46-47: Find out how many plot panels are needed; implemented by `if ismember('spectrogram',spin_system.control.plotting)`.
- Lines 60-61: Get a figure without stealing focus; implemented by `try`.
- Lines 68-69: Plot spectrograms for as long as timing grid permits; implemented by `if ismember('spectrogram',spin_system.control.plotting)`.
- Lines 71-72: Count uniform slices; implemented by `last_slice_to_plot=1`.
- Lines 82-83: Refuse when too few; implemented by `if last_slice_to_plot<5`.
- Lines 87-88: Truncate the timing grid for spectrogram plotting; implemented by `timing_grid=spin_system.control.pulse_dt(1:last_slice_to_plot)`.
- Lines 90-91: Get the sampling rate; implemented by `sampl_rate=1/spin_system.control.pulse_dt(1)`.
- Lines 93-94: Loop over channel pairs; implemented by `for n=1:(size(waveform,1)/2)`.
- Lines 96-97: Set the current plot; implemented by `subplot(n_plots_x,n_plots_y,current_plot)`.
- Lines 99-101: Get the complex channel waveform; implemented by `cplx_ch_wf= waveform(2*n-1,1:last_slice_to_plot) -1i*waveform(2*n ,1:last_slice_to_plot)`.
- Lines 103-104: Mirror replication at the edges; implemented by `padded_wf=[fliplr(cplx_ch_wf) cplx_ch_wf fliplr(cplx_ch_wf)]`.
- Lines 106-107: Optimal spectrogram settings; implemented by `n_slices=numel(timing_grid)`.
- Lines 113-115: Get the spectrogram at the optimal visually informative settings; implemented by `[st_fft,f_axis,t_axis]=spectrogram(padded_wf,window_size,window_overlap, n_freq_bins,sampl_rate,'yaxis','center')`.

### Control flow inferred from the code

- Line 32: conditional branch on `isempty(spin_system.control.plotting), return; end`.
- Line 44: conditional branch on `n_plots==0, return; end`.
- Line 47: conditional branch on `ismember('spectrogram',spin_system.control.plotting)`.
- Line 50: conditional branch on `ismember('time_by_slice',spin_system.control.plotting)`.
- Line 53: conditional branch on `n_plots>3`.
- Line 69: conditional branch on `ismember('spectrogram',spin_system.control.plotting)`.
- Line 73: `for` loop over `n=1:numel(spin_system.control.pulse_dt)`.
- Line 74: conditional branch on `spin_system.control.pulse_dt(n)==`.
- Line 83: conditional branch on `last_slice_to_plot<5`.
- Line 94: `for` loop over `n=1:(size(waveform,1)/2)`.
- Line 136: conditional branch on `last_slice_to_plot==numel(spin_system.control.pulse_dt)`.
- Line 156: conditional branch on `ismember('time_by_slice',spin_system.control.plotting)`.
- Line 174: dispatches on `spin_system.control.integrator`; cases `'rectangle'`, `'trapezium'`, `'trapezium'`, `'rectangle'`.
- Line 205: conditional branch on `strcmp(spin_system.control.integrator,'rectangle')&&(~isscalar(lower_bound))`.

### Key state/data transformations

- Lines 38: computes `spin_system.control.plotting` using `spin_system.control.plotting=setdiff(spin_system.control.plotting,{'trajectory'})`.
- Lines 41: computes `n_plots` using `n_plots=numel(spin_system.control.plotting)`.
- Lines 54: computes `n_plots_y` using `n_plots_y=ceil(n_plots/2); n_plots_x=2`.
- Lines 58: computes `current_plot` using `current_plot=1`.
- Lines 72: computes `last_slice_to_plot` using `last_slice_to_plot=1`.
- Lines 88: computes `timing_grid` using `timing_grid=spin_system.control.pulse_dt(1:last_slice_to_plot)`.
- Lines 91: computes `sampl_rate` using `sampl_rate=1/spin_system.control.pulse_dt(1)`.
- Lines 100-101: computes `cplx_ch_wf` using `cplx_ch_wf= waveform(2*n-1,1:last_slice_to_plot) -1i*waveform(2*n ,1:last_slice_to_plot)`.
- Lines 104: computes `padded_wf` using `padded_wf=[fliplr(cplx_ch_wf) cplx_ch_wf fliplr(cplx_ch_wf)]`.
- Lines 107: computes `n_slices` using `n_slices=numel(timing_grid)`.
- Lines 108: computes `window_size` using `window_size=min(n_slices,ceil(sqrt(2*n_slices)))`.
- Lines 109: computes `window_step` using `window_step=max(1,ceil(window_size/4))`.
- Lines 110: computes `window_overlap` using `window_overlap=window_size-window_step`.
- Lines 111: computes `n_freq_bins` using `n_freq_bins=2*window_size`.
- Lines 114-115: computes `[st_fft,f_axis,t_axis]` using `[st_fft,f_axis,t_axis]=spectrogram(padded_wf,window_size,window_overlap, n_freq_bins,sampl_rate,'yaxis','center')`.
- Lines 118: computes `t_axis` using `t_axis=t_axis-sum(timing_grid)`.
- Lines 121: computes `phi` using `phi=atan2(real(st_fft),imag(st_fft)); phi=(phi+pi)/(2*pi)`.
- Lines 122: computes `amp` using `amp=abs(st_fft); amp=amp/max(amp,[],'all')`.

### Local helper functions

- Line 755: `grumble()` — `function grumble(spin_system,waveform,traj_data,fidelities)`.
  - Representative operation: `if (~isnumeric(waveform))||(~isreal(waveform))`.
  - Representative operation: `error('waveform must be a real numerical array.')`.

## Parameters / inputs

- waveform -waveform, as supplied to a user-end
- function, such as grape_xy.m
- traj_data -a cell array of trajectory data struc-
- tures returned by GRAPE, one per en-
- semble member
- fidelities -fidelities array, as returned by
- user-end functions, such as grape_xy
- Note: this function is called internally by the optimal cont-
- rol module, you should not be calling it directly. All
- settings should be specified in the call to optimcon.m
- when the optimal control problem is set up.

## Implementation structure

- Diagnostic plotting function for optimal control module. Plots
- trajectory and control pulse analysis. Syntax:
- ctrl_trajan(spin_system,waveform,trajectory,fidelities)
- waveform -waveform, as supplied to a user-end
- function, such as grape_xy.m
- traj_data -a cell array of trajectory data struc-
- tures returned by GRAPE, one per en-
- semble member
- fidelities -fidelities array, as returned by
- user-end functions, such as grape_xy
- Note: this function is called internally by the optimal cont-
- rol module, you should not be calling it directly. All

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `setdiff()`, `ismember()`, `set()`, `scale_figure()`, `subplot()`, `waveform()`, `fliplr()`, `spectrogram()`, `atan2()`, `cat()`, `image()`, `hsv2rgb()`, `ktitle()`, `num2str()`, `kylabel()`.
