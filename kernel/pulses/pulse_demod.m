% Interactive demodulation of a complex pulse waveform by a user-
% specified frequency. Syntax:
%
%              demod_pulse=pulse_demod(time_grid,in_phase,out_phase)
%
% Parameters:
%
%   time_grid  - strictly increasing time grid, seconds
%
%   in_phase   - in-phase pulse waveform component
%
%   out_phase  - out-of-phase pulse waveform component
%
% Outputs:
%
%   demod_pulse - demodulated complex pulse waveform
%
% Notes:
%
%   The frequency entry field uses Hz. The GHz, MHz, kHz, and Hz buttons
%   set the slider step size without changing the demodulation frequency.
%   The phase and frequency buttons switch the plot between unwrapped
%   phase in radians and instantaneous frequency in Hz. The sticky wrap
%   button switches the phase plot into the [0,2*pi] interval.
%   The slider range is a moving 100-step window in the selected units.
%   The save button returns the current demodulated waveform and exits.
%   The complex waveform in_phase+1i*out_phase is multiplied by
%   exp(2*pi*1i*freq*time_grid).
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pulse_demod.m>

function demod_pulse=pulse_demod(time_grid,in_phase,out_phase)

% Check consistency
grumble(time_grid,in_phase,out_phase);

% Build the complex pulse waveform
complex_pulse=in_phase+1i*out_phase;

% Initialise the returned waveform
demod_pulse=complex_pulse;

% Start at zero-frequency demodulation
freq=0; freq_step=1; step_name='Hz';
plot_mode='phase'; phase_wrapped=false;
shown_mode='phase'; shown_wrapped=false;

% Set the initial slider range from the active step
range_steps=100;
freq_low=-range_steps*freq_step/2;
freq_high=+range_steps*freq_step/2;

% Create the figure window
fig_handle=kfigure('Name','Pulse demodulation',...
                   'NumberTitle','off');
scale_figure([1.4 1.2]);

% Create the waveform axis
axis_handle=axes('Parent',fig_handle,...
                 'Units','normalized',...
                 'Position',[0.10 0.40 0.70 0.54]);

% Plot the initial waveform representation
plot_line=plot(axis_handle,time_grid,zeros(size(time_grid)),'b-');
axis(axis_handle,'tight'); kgrid;
ktitle('demodulation frequency: 0 Hz');
kxlabel('time, s'); kylabel('phase, rad');
legend_handle=klegend(axis_handle,plot_line,...
                      {'unwrapped phase'},...
                      'Location','Best');
title_handle=get(axis_handle,'Title');
ylabel_handle=get(axis_handle,'YLabel');

% Create the right-side display buttons
button_phase=uicontrol('Parent',fig_handle,'Style','pushbutton',...
                       'Units','normalized',...
                       'Position',[0.84 0.82 0.12 0.06],...
                       'String','phase',...
                       'Callback',{@set_plot,'phase'});
button_freq=uicontrol('Parent',fig_handle,'Style','pushbutton',...
                      'Units','normalized',...
                      'Position',[0.84 0.74 0.12 0.06],...
                      'String','frequency',...
                      'Callback',{@set_plot,'frequency'});
button_wrap=uicontrol('Parent',fig_handle,'Style','togglebutton',...
                      'Units','normalized',...
                      'Position',[0.84 0.66 0.12 0.06],...
                      'String','wrap',...
                      'Callback',@set_wrap);

% Label the slider step controls
uicontrol('Parent',fig_handle,'Style','text',...
          'Units','normalized',...
          'Position',[0.28 0.27 0.54 0.035],...
          'String','slider step',...
          'HorizontalAlignment','left');

% Create the order-of-magnitude buttons
button_ghz=uicontrol('Parent',fig_handle,'Style','pushbutton',...
                     'Units','normalized',...
                     'Position',[0.28 0.20 0.12 0.06],...
                     'String','GHz',...
                     'Callback',{@set_mult,1e9,'GHz'});
button_mhz=uicontrol('Parent',fig_handle,'Style','pushbutton',...
                     'Units','normalized',...
                     'Position',[0.42 0.20 0.12 0.06],...
                     'String','MHz',...
                     'Callback',{@set_mult,1e6,'MHz'});
button_khz=uicontrol('Parent',fig_handle,'Style','pushbutton',...
                     'Units','normalized',...
                     'Position',[0.56 0.20 0.12 0.06],...
                     'String','kHz',...
                     'Callback',{@set_mult,1e3,'kHz'});
button_hz=uicontrol('Parent',fig_handle,'Style','pushbutton',...
                    'Units','normalized',...
                    'Position',[0.70 0.20 0.12 0.06],...
                    'String','Hz',...
                    'Callback',{@set_mult,1,'Hz'});

% Create the save button
uicontrol('Parent',fig_handle,'Style','pushbutton',...
          'Units','normalized',...
          'Position',[0.84 0.20 0.12 0.06],...
          'String','save',...
          'Callback',@save_pulse);

% Create the frequency entry field
edit_handle=uicontrol('Parent',fig_handle,'Style','edit',...
                      'Units','normalized',...
                      'Position',[0.10 0.10 0.14 0.06],...
                      'String','0',...
                      'Callback',@edit_freq);

% Label the frequency entry field
uicontrol('Parent',fig_handle,'Style','text',...
          'Units','normalized',...
          'Position',[0.10 0.16 0.14 0.035],...
          'String','frequency, Hz',...
          'HorizontalAlignment','left');

% Set the slider step from the initial multiplier
minor_step=max(min(freq_step/(freq_high-freq_low),1),eps);
major_step=max(min(10*minor_step,1),minor_step);

% Create the frequency slider
slider_handle=uicontrol('Parent',fig_handle,'Style','slider',...
                        'Units','normalized',...
                        'Position',[0.28 0.10 0.68 0.06],...
                        'Min',freq_low,...
                        'Max',freq_high,...
                        'Value',0,...
                        'SliderStep',[minor_step major_step],...
                        'Callback',@slide_freq);

% Highlight the initial slider step and display mode
select_button(); select_plot(); select_wrap();
update_plot();

% Wait until the user saves or closes the figure
uiwait(fig_handle);

% Close the figure after a save request
if ishandle(fig_handle)
    delete(fig_handle);
end

    % Change the active slider step size
    function set_mult(~,~,new_mult,new_name)

        % Store the new slider step size
        freq_step=new_mult; step_name=new_name;

        % Update the slider without changing the frequency
        set_slider(true);

        % Update the active button and plot
        select_button(); update_plot();

    end

    % Change the active plot mode
    function set_plot(~,~,new_mode)

        % Store the new plot mode
        plot_mode=new_mode;

        % Update the active button and plot
        select_plot(); update_plot();

    end

    % Change the active phase wrapping state
    function set_wrap(~,~)

        % Store the wrapping button state
        phase_wrapped=get(button_wrap,'Value')>0;

        % Update the button text and plot
        select_wrap(); update_plot();

    end

    % Read frequency changes from the slider
    function slide_freq(~,~)

        % Update the demodulation frequency
        freq=get(slider_handle,'Value');

        % Update the frequency display
        set(edit_handle,'String',num2str(freq,'%.12g'));
        set_slider(false); update_plot();

    end

    % Read frequency changes from the entry field
    function edit_freq(~,~)

        % Interpret the entry in Hz
        edit_value=str2double(get(edit_handle,'String'));

        % Reject invalid frequencies
        if ~isfinite(edit_value)
            set(edit_handle,'String',num2str(freq,'%.12g'));
            return
        end

        % Update the demodulation frequency
        freq=edit_value;

        % Update the slider and plot
        set_slider(true); update_plot();

    end

    % Update the slider range and step size
    function set_slider(centre_range)

        % Set the slider window size in active units
        target_span=range_steps*freq_step;

        % Measure the current slider range
        freq_span=freq_high-freq_low;

        % Recentre the window after typed or unit changes
        if centre_range||(freq_span>target_span)
            freq_low=freq-target_span/2;
            freq_high=freq+target_span/2;
            freq_span=target_span;
        end

        % Shift the lower end when the thumb approaches it
        if (freq-freq_low)<10*freq_step
            freq_low=freq-10*freq_step;
            freq_high=freq_low+freq_span;
        end

        % Shift the upper end when the thumb approaches it
        if (freq_high-freq)<10*freq_step
            freq_high=freq+10*freq_step;
            freq_low=freq_high-freq_span;
        end

        % Set one small slider step to one multiplier unit
        minor_step=max(min(freq_step/(freq_high-freq_low),1),eps);
        major_step=max(min(10*minor_step,1),minor_step);

        % Update slider geometry in Hz
        set(slider_handle,'Min',freq_low,...
                          'Max',freq_high,...
                          'Value',freq,...
                          'SliderStep',[minor_step major_step]);

    end

    % Show the active slider step size
    function select_button()

        % Clear all button highlights
        set([button_ghz button_mhz button_khz button_hz],...
            'FontWeight','normal');

        % Highlight the active button
        switch step_name

            case 'GHz'

                set(button_ghz,'FontWeight','bold');

            case 'MHz'

                set(button_mhz,'FontWeight','bold');

            case 'kHz'

                set(button_khz,'FontWeight','bold');

            case 'Hz'

                set(button_hz,'FontWeight','bold');

        end

    end

    % Show the active plot mode
    function select_plot()

        % Clear all button highlights
        set([button_phase button_freq],'FontWeight','normal');

        % Highlight the active button
        switch plot_mode

            case 'phase'

                set(button_phase,'FontWeight','bold');

            case 'frequency'

                set(button_freq,'FontWeight','bold');

        end

    end

    % Show the active phase wrapping state
    function select_wrap()

        % Label the button with the available action
        if phase_wrapped
            set(button_wrap,'String','unwrap',...
                            'Value',1,...
                            'FontWeight','bold');
        else
            set(button_wrap,'String','wrap',...
                            'Value',0,...
                            'FontWeight','normal');
        end

    end

    % Save the current waveform and leave the GUI
    function save_pulse(~,~)

        % Accept any pending frequency-field edit
        edit_freq([],[]);

        % Leave the modal wait state
        uiresume(fig_handle);

    end

    % Redraw the demodulated waveform
    function update_plot()

        % Apply the demodulation multiplier
        demod_pulse=complex_pulse.*exp(2*pi*1i*freq*time_grid);

        % Compute the demodulated phase
        pulse_phase=unwrap(angle(demod_pulse));

        % Display the selected waveform representation
        switch plot_mode

            case 'phase'

                % Apply phase wrapping if requested
                if phase_wrapped
                    plot_phase=mod(pulse_phase,2*pi);
                    legend_text='wrapped phase';
                else
                    plot_phase=pulse_phase;
                    legend_text='unwrapped phase';
                end

                set(plot_line,'XData',time_grid,...
                              'YData',plot_phase);
                if (~strcmp(shown_mode,'phase'))||...
                   (shown_wrapped~=phase_wrapped)
                    set(ylabel_handle,'String','phase, rad');
                    delete(legend_handle);
                    legend_handle=klegend(axis_handle,plot_line,...
                                          {legend_text},...
                                          'Location','Best');
                    shown_mode='phase'; shown_wrapped=phase_wrapped;
                end

            case 'frequency'

                freq_grid=(time_grid(1:(end-1))+time_grid(2:end))/2;
                inst_frq=diff(pulse_phase)./(2*pi*diff(time_grid));
                set(plot_line,'XData',freq_grid,...
                              'YData',inst_frq);
                if ~strcmp(shown_mode,'frequency')
                    set(ylabel_handle,'String','instantaneous frequency, Hz');
                    delete(legend_handle);
                    legend_handle=klegend(axis_handle,plot_line,...
                                          {'instantaneous frequency'},...
                                          'Location','Best');
                    shown_mode='frequency';
                end

        end

        % Update the plot annotations
        set(title_handle,'String',['\textbf{demodulation frequency: '...
                         num2str(freq,'%.12g') ' Hz}']);
        axis(axis_handle,'tight'); drawnow nocallbacks;

    end

end

% Consistency enforcement
function grumble(time_grid,in_phase,out_phase)
if (~isnumeric(time_grid))||(~isreal(time_grid))||(~isvector(time_grid))||...
   (numel(time_grid)<2)||any(~isfinite(time_grid))
    error('time_grid must be a finite real vector with at least two elements.');
end
if any(diff(time_grid)<=0)
    error('time_grid must be strictly increasing.');
end
if (~isnumeric(in_phase))||(~isreal(in_phase))||(~isvector(in_phase))||...
   any(~isfinite(in_phase))
    error('in_phase must be a finite real vector.');
end
if (~isnumeric(out_phase))||(~isreal(out_phase))||(~isvector(out_phase))||...
   any(~isfinite(out_phase))
    error('out_phase must be a finite real vector.');
end
if (~all(size(time_grid)==size(in_phase)))||...
   (~all(size(in_phase)==size(out_phase)))
    error('time_grid, in_phase, and out_phase must have the same dimension.');
end
end

% Gods do not like happy people.
%
% Herodotus

