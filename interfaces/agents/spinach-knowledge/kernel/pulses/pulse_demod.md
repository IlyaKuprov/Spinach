# kernel/pulses/pulse_demod.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/pulse_demod.m`
- Signature: `demod_pulse=pulse_demod(time_grid,in_phase,out_phase)`
- Total lines: 442

## Purpose

Interactive demodulation of a complex pulse waveform by a user- specified frequency. Syntax: demod_pulse=pulse_demod(time_grid,in_phase,out_phase)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `set_mult()`, `set_plot()`, `get()`, `str2double()`, `set_slider()`, `select_button()`, `select_plot()`, `select_wrap()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(time_grid,in_phase,out_phase)`.
- Lines 39-40: Build the complex pulse waveform; implemented by `complex_pulse=in_phase+1i*out_phase`.
- Lines 42-43: Initialise the returned waveform; implemented by `demod_pulse=complex_pulse`.
- Lines 45-46: Start at zero-frequency demodulation; implemented by `freq=0; freq_step=1; step_name='Hz'`.
- Lines 50-51: Set the initial slider range from the active step; implemented by `range_steps=100`.
- Lines 55-57: Create the figure window; implemented by `fig_handle=kfigure('Name','Pulse demodulation', 'NumberTitle','off')`.
- Lines 60-63: Create the waveform axis; implemented by `axis_handle=axes('Parent',fig_handle, 'Units','normalized', 'Position',[0.10 0.40 0.70 0.54])`.
- Lines 65-66: Plot the initial waveform representation; implemented by `plot_line=plot(axis_handle,time_grid,zeros(size(time_grid)),'b-')`.
- Lines 76-81: Create the right-side display buttons; implemented by `button_phase=uicontrol('Parent',fig_handle,'Style','pushbutton', 'Units','normalized', 'Position',[0.84 0.82 0.12 0.06], 'String','phase', 'Callback',{@set_plot,'phase'})`.
- Lines 93-98: Label the slider step controls; implemented by `uicontrol('Parent',fig_handle,'Style','text', 'Units','normalized', 'Position',[0.28 0.27 0.54 0.035], 'String','slider step', 'HorizontalAlignment','left')`.
- Lines 100-105: Create the order-of-magnitude buttons; implemented by `button_ghz=uicontrol('Parent',fig_handle,'Style','pushbutton', 'Units','normalized', 'Position',[0.28 0.20 0.12 0.06], 'String','GHz', 'Callback',{@set_mult,1e9,'GHz'})`.
- Lines 122-127: Create the save button; implemented by `uicontrol('Parent',fig_handle,'Style','pushbutton', 'Units','normalized', 'Position',[0.84 0.20 0.12 0.06], 'String','save', 'Callback',@save_pulse)`.
- Lines 129-134: Create the frequency entry field; implemented by `edit_handle=uicontrol('Parent',fig_handle,'Style','edit', 'Units','normalized', 'Position',[0.10 0.10 0.14 0.06], 'String','0', 'Callback',@edit_freq)`.
- Lines 136-141: Label the frequency entry field; implemented by `uicontrol('Parent',fig_handle,'Style','text', 'Units','normalized', 'Position',[0.10 0.16 0.14 0.035], 'String','frequency, Hz', 'HorizontalAlignment','left')`.
- Lines 143-144: Set the slider step from the initial multiplier; implemented by `minor_step=max(min(freq_step/(freq_high-freq_low),1),eps)`.
- Lines 147-155: Create the frequency slider; implemented by `slider_handle=uicontrol('Parent',fig_handle,'Style','slider', 'Units','normalized', 'Position',[0.28 0.10 0.68 0.06], 'Min',freq_low, 'Max',freq_high, 'Value',0, 'Slider…`.
- Lines 157-158: Highlight the initial slider step and display mode; implemented by `select_button(); select_plot(); select_wrap()`.
- Lines 161-162: Wait until the user saves or closes the figure; implemented by `uiwait(fig_handle)`.

### Control flow inferred from the code

- Line 165: conditional branch on `ishandle(fig_handle)`.

### Key state/data transformations

- Lines 40: computes `complex_pulse` using `complex_pulse=in_phase+1i*out_phase`.
- Lines 43: computes `demod_pulse` using `demod_pulse=complex_pulse`.
- Lines 46: computes `freq` using `freq=0; freq_step=1; step_name='Hz'`.
- Lines 47: computes `plot_mode` using `plot_mode='phase'; phase_wrapped=false`.
- Lines 48: computes `shown_mode` using `shown_mode='phase'; shown_wrapped=false`.
- Lines 51: computes `range_steps` using `range_steps=100`.
- Lines 52: computes `freq_low` using `freq_low=-range_steps*freq_step/2`.
- Lines 53: computes `freq_high` using `freq_high=+range_steps*freq_step/2`.
- Lines 56-57: computes `fig_handle` using `fig_handle=kfigure('Name','Pulse demodulation', 'NumberTitle','off')`.
- Lines 61-63: computes `axis_handle` using `axis_handle=axes('Parent',fig_handle, 'Units','normalized', 'Position',[0.10 0.40 0.70 0.54])`.
- Lines 66: computes `plot_line` using `plot_line=plot(axis_handle,time_grid,zeros(size(time_grid)),'b-')`.
- Lines 70-72: computes `legend_handle` using `legend_handle=klegend(axis_handle,plot_line, {'unwrapped phase'}, 'Location','Best')`.
- Lines 73: computes `title_handle` using `title_handle=get(axis_handle,'Title')`.
- Lines 74: computes `ylabel_handle` using `ylabel_handle=get(axis_handle,'YLabel')`.
- Lines 77-81: computes `button_phase` using `button_phase=uicontrol('Parent',fig_handle,'Style','pushbutton', 'Units','normalized', 'Position',[0.84 0.82 0.12 0.06], 'String','phase', 'Callback',{@set_plot,'phase'})`.
- Lines 82-86: computes `button_freq` using `button_freq=uicontrol('Parent',fig_handle,'Style','pushbutton', 'Units','normalized', 'Position',[0.84 0.74 0.12 0.06], 'String','frequency', 'Callback',{@set_plot,'freq…`.
- Lines 87-91: computes `button_wrap` using `button_wrap=uicontrol('Parent',fig_handle,'Style','togglebutton', 'Units','normalized', 'Position',[0.84 0.66 0.12 0.06], 'String','wrap', 'Callback',@set_wrap)`.
- Lines 101-105: computes `button_ghz` using `button_ghz=uicontrol('Parent',fig_handle,'Style','pushbutton', 'Units','normalized', 'Position',[0.28 0.20 0.12 0.06], 'String','GHz', 'Callback',{@set_mult,1e9,'GHz'})`.

### Local helper functions

- Line 170: `set_mult()` — `function set_mult(~,~,new_mult,new_name)`. Store the new slider step size
  - Representative operation: `freq_step=new_mult; step_name=new_name`.
  - Representative operation: `set_slider(true)`.
- Line 184: `set_plot()` — `function set_plot(~,~,new_mode)`. Store the new plot mode
  - Representative operation: `plot_mode=new_mode`.
  - Representative operation: `select_plot(); update_plot()`.
- Line 195: `set_wrap()` — `function set_wrap(~,~)`. Store the wrapping button state
  - Representative operation: `phase_wrapped=get(button_wrap,'Value')>0`.
  - Representative operation: `select_wrap(); update_plot()`.
- Line 206: `slide_freq()` — `function slide_freq(~,~)`. Update the demodulation frequency
  - Representative operation: `freq=get(slider_handle,'Value')`.
  - Representative operation: `set(edit_handle,'String',num2str(freq,'%.12g'))`.
- Line 218: `edit_freq()` — `function edit_freq(~,~)`. Interpret the entry in Hz
  - Representative operation: `edit_value=str2double(get(edit_handle,'String'))`.
  - Representative operation: `if ~isfinite(edit_value)`.
- Line 238: `set_slider()` — `function set_slider(centre_range)`. Set the slider window size in active units
  - Representative operation: `target_span=range_steps*freq_step`.
  - Representative operation: `freq_span=freq_high-freq_low`.
- Line 278: `select_button()` — `function select_button()`. Clear all button highlights
  - Representative operation: `set([button_ghz button_mhz button_khz button_hz], 'FontWeight','normal')`.
  - Representative operation: `'FontWeight','normal')`.
- Line 308: `select_plot()` — `function select_plot()`. Clear all button highlights
  - Representative operation: `set([button_phase button_freq],'FontWeight','normal')`.
  - Representative operation: `switch plot_mode`.

## Parameters / inputs

- time_grid -strictly increasing time grid, seconds
- in_phase -in-phase pulse waveform component
- out_phase -out-of-phase pulse waveform component

## Outputs

- demod_pulse -demodulated complex pulse waveform

## Header notes

- The frequency entry field uses Hz. The GHz, MHz, kHz, and Hz buttons
- set the slider step size without changing the demodulation frequency.
- The phase and frequency buttons switch the plot between unwrapped
- phase in radians and instantaneous frequency in Hz. The sticky wrap
- button switches the phase plot into the [0,2*pi] interval.
- The slider range is a moving 100-step window in the selected units.
- The save button returns the current demodulated waveform and exits.
- The complex waveform in_phase+1i*out_phase is multiplied by
- exp(2*pi*1i*freq*time_grid).

## Implementation structure

- Interactive demodulation of a complex pulse waveform by a user-
- specified frequency. Syntax:
- demod_pulse=pulse_demod(time_grid,in_phase,out_phase)
- time_grid -strictly increasing time grid, seconds
- in_phase -in-phase pulse waveform component
- out_phase -out-of-phase pulse waveform component
- demod_pulse -demodulated complex pulse waveform
- The frequency entry field uses Hz. The GHz, MHz, kHz, and Hz buttons
- set the slider step size without changing the demodulation frequency.
- The phase and frequency buttons switch the plot between unwrapped
- phase in radians and instantaneous frequency in Hz. The sticky wrap
- button switches the phase plot into the [0,2*pi] interval.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `kfigure()`, `scale_figure()`, `axes()`, `ktitle()`, `kxlabel()`, `kylabel()`, `klegend()`, `get()`, `uicontrol()`, `select_button()`, `select_plot()`, `select_wrap()`, `update_plot()`, `uiwait()`, `ishandle()`.
