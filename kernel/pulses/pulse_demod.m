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
%   set the slider multiplier without changing the demodulation frequency.
%   The slider range expands when its end stops are reached.
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

% Set the initial slider range from the sampling limit
init_range=1/(2*min(diff(time_grid)));
freq_low=-init_range; freq_high=+init_range;

% Start at zero-frequency demodulation
freq=0; slider_mult=1; slider_name='Hz';

% Create the figure window
fig_handle=kfigure('Name','Pulse demodulation',...
                   'NumberTitle','off');
scale_figure([1.4 1.2]);

% Create the waveform axis
axis_handle=axes('Parent',fig_handle,...
                 'Units','normalized',...
                 'Position',[0.10 0.40 0.86 0.54]);

% Plot the initial waveform components
in_line=plot(axis_handle,time_grid,in_phase,'b-');
hold(axis_handle,'on');
out_line=plot(axis_handle,time_grid,out_phase,'r-');
axis(axis_handle,'tight'); kgrid;
ktitle('demodulation frequency: 0 Hz');
kxlabel('time, s'); kylabel('pulse amplitude');
klegend(axis_handle,{'in-phase','out-of-phase'},'Location','best');
title_handle=get(axis_handle,'Title');

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

% Set the slider step from the initial multiplier
minor_step=max(min(slider_mult/(freq_high-freq_low),1),eps);
major_step=max(min(10*minor_step,1),minor_step);

% Create the frequency slider
slider_handle=uicontrol('Parent',fig_handle,'Style','slider',...
                        'Units','normalized',...
                        'Position',[0.28 0.10 0.68 0.06],...
                        'Min',freq_low/slider_mult,...
                        'Max',freq_high/slider_mult,...
                        'Value',0,...
                        'SliderStep',[minor_step major_step],...
                        'Callback',@slide_freq);

% Highlight the initial slider multiplier
select_button();

% Wait until the user saves or closes the figure
uiwait(fig_handle);

% Close the figure after a save request
if ishandle(fig_handle)
    delete(fig_handle);
end


    % Change the active slider multiplier
    function set_mult(~,~,new_mult,new_name)

        % Store the new multiplier
        slider_mult=new_mult; slider_name=new_name;

        % Rescale the slider without changing the frequency
        set_slider();

        % Update the active button and plot
        select_button(); update_plot();

    end

    % Read frequency changes from the slider
    function slide_freq(~,~)

        % Update the demodulation frequency
        freq=get(slider_handle,'Value')*slider_mult;

        % Measure the current slider range
        freq_span=freq_high-freq_low;

        % Create upper headroom if the upper stop is reached
        if freq>=freq_high
            freq_high=freq+freq_span;
        end

        % Create lower headroom if the lower stop is reached
        if freq<=freq_low
            freq_low=freq-freq_span;
        end

        % Update the frequency display
        set(edit_handle,'String',num2str(freq,'%.12g'));
        set_slider(); update_plot();

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

        % Measure the current slider range
        freq_span=freq_high-freq_low;

        % Create upper headroom if the upper stop is reached
        if freq>=freq_high
            freq_high=freq+freq_span;
        end

        % Create lower headroom if the lower stop is reached
        if freq<=freq_low
            freq_low=freq-freq_span;
        end

        % Update the slider and plot
        set_slider(); update_plot();

    end

    % Update the slider range and step size
    function set_slider()

        % Set one small slider step to one multiplier unit
        minor_step=max(min(slider_mult/(freq_high-freq_low),1),eps);
        major_step=max(min(10*minor_step,1),minor_step);

        % Update slider geometry without changing the frequency
        set(slider_handle,'Min',freq_low/slider_mult,...
                          'Max',freq_high/slider_mult,...
                          'Value',freq/slider_mult,...
                          'SliderStep',[minor_step major_step]);

    end

    % Show the active slider multiplier
    function select_button()

        % Clear all button highlights
        set([button_ghz button_mhz button_khz button_hz],...
            'FontWeight','normal');

        % Highlight the active button
        switch slider_name

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

        % Update the in-phase and out-of-phase traces
        set(in_line,'YData',real(demod_pulse));
        set(out_line,'YData',imag(demod_pulse));

        % Update the plot annotations
        set(title_handle,'String',['\textbf{demodulation frequency: '...
                         num2str(freq,'%.12g') ' Hz}']);
        axis(axis_handle,'tight'); drawnow;

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


