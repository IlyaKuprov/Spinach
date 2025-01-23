% Diagnostic plotting function for optimal control module. Plots
% trajectory and control pulse analysis. Syntax:
%
%     ctrl_trajan(spin_system,waveform,trajectory,fidelities)
%
% Parameters:
%
%    waveform   - waveform, as supplied to a user-end
%                 function, such as grape_xy.m
%
%    traj_data  - a cell array of trajectory data struc-
%                 tures returned by GRAPE, one per en-
%                 semble member
%
%    fidelities - fidelities array, as returned by
%                 user-end functions, such as grape_xy
%
% Note: this function is called internally by the optimal cont-
%       rol module, you should not be calling it directly. All
%       settings should be specified in the call to optimcon.m
%       when the optimal control problem is set up.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ctrl_trajan.m>

function ctrl_trajan(spin_system,waveform,traj_data,fidelities)

% If no plotting is needed, exit
if isempty(spin_system.control.plotting), return; end

% Check consistency
grumble(spin_system,waveform,traj_data,fidelities);

% Count the plots
n_plots=numel(spin_system.control.plotting); 
if ismember('spectrogram',spin_system.control.plotting)
    n_plots=n_plots+size(waveform,1)/2-1;
end
if n_plots>3
    n_plots_y=ceil(n_plots/2); n_plots_x=2;
else
    n_plots_y=n_plots; n_plots_x=1;
end
current_plot=1;

% Get a figure without stealing focus
try
    set(groot,'CurrentFigure',1);
    scale_figure([0.70*n_plots_y 0.70*n_plots_x]);
catch
    figure(1); set(groot,'CurrentFigure',1);
    scale_figure([0.70*n_plots_y 0.70*n_plots_x]);
end

% Plot the spectrograms
if ismember('spectrogram',spin_system.control.plotting)

    % Loop over channel pairs
    for n=1:(size(waveform,1)/2)

        % Set the current plot
        subplot(n_plots_x,n_plots_y,current_plot);

        % Get the complex channel waveform
        cplx_ch_wf=waveform(2*n-1,:)-1i*waveform(2*n,:);

        % Mirror replication at the edges
        padded_wf=[fliplr(cplx_ch_wf) cplx_ch_wf fliplr(cplx_ch_wf)];

        % Get the sampling rate
        sampl_rate=1/spin_system.control.pulse_dt(1);

        % Get the spectrogram
        window_size=ceil(sqrt(numel(spin_system.control.pulse_dt)));
        window_overlap=ceil(window_size/2); n_freq_bins=2*window_size;
        [st_fft,f_axis,t_axis]=spectrogram(padded_wf,window_size,window_overlap,...
                                           n_freq_bins,sampl_rate,'yaxis','center');

        % Modify the time axis to account for the replicas
        t_axis=t_axis-sum(spin_system.control.pulse_dt);
        
        % Interpret phase as hue and amplitude as value in HSV - phase unwrap needed
        phi=atan2(real(st_fft),imag(st_fft)); phi=(phi+pi)/(2*pi); 
        amp=abs(st_fft); amp=amp/max(amp,[],'all');
        hsv=cat(3,ones(size(phi)),ones(size(phi)),amp); rgb=hsv2rgb(hsv);
        image(t_axis,f_axis,rgb); kxlabel('time, seconds');
        ktitle(['channels ' num2str(2*n-1) ',' num2str(2*n)]);
        kylabel('frequency offset, Hz'); set(gca,'YDir','normal');

        % Only keep the physically relevan interval
        xlim([0 sum(spin_system.control.pulse_dt)]);

        % Increment plot number
        current_plot=current_plot+1;

    end

end

% Build the time axis
t_axis=[0 cumsum(spin_system.control.pulse_dt)];

% Apply the average power for plotting
average_power=mean(spin_system.control.pwr_levels);

% Adapt to the integrator type
switch spin_system.control.integrator
    
    % Piecewise-constant
    case 'rectangle'

        % Compute maximum nutation angle in degrees
        max_nut=max(abs(spin_system.control.pulse_dt.*waveform),[],'all');
        max_nut=180*average_power*max_nut/pi;
        
        % An extra point at the end for stairs plot
        waveform=average_power*[waveform waveform(:,end)];

    % Piecewise-linear
    case 'trapezium'
        
        % Compute maximum nutation angle in degrees
        medians=(waveform(:,2:end)+...
                 waveform(:,1:(end-1)))/2;
        max_nut=max(abs(spin_system.control.pulse_dt.*medians),[],'all');
        max_nut=180*average_power*max_nut/pi;

        % Just apply the power level
        waveform=average_power*waveform;

end

% Compute upper and lower power bounds
lower_bound=spin_system.control.l_bound*average_power;
upper_bound=spin_system.control.u_bound*average_power;

% Plot Cartesian controls
if ismember('xy_controls',spin_system.control.plotting)
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot);
    
    % Adapt to the integrator type
    switch spin_system.control.integrator
        
        % Piecewise-linear
        case 'trapezium'

            % Conventional linear plot
            p=plot(t_axis',waveform'/(2*pi));

        % Piecewise-constant
        case 'rectangle'
            
            % Stairs plot
            p=stairs(t_axis',waveform'/(2*pi));

        otherwise

            % Complain and bomb out
            error('unknown integrator.');

    end

    % Plot with predictable colours
    for n=1:numel(p)
        p(n).Color=hsv2rgb([n/numel(p) 0.75 0.75]);
    end
    
    % Power bounds
    fb_pwr=refline([0 lower_bound/(2*pi)]);
    set(fb_pwr,'Color','k','LineStyle','--');
    cb_pwr=refline([0 upper_bound/(2*pi)]);
    set(cb_pwr,'Color','k','LineStyle','--');
    
    % Compute axis limits
    x_lower=0; x_upper=spin_system.control.pulse_dur;
    y_lower=min([min(waveform(:)), lower_bound-0.05*abs(lower_bound)])/(2*pi);
    y_upper=max([max(waveform(:)), upper_bound+0.05*abs(upper_bound)])/(2*pi);

    % Set axis limits
    xlim([x_lower x_upper]); ylim([y_lower y_upper]);
       
    % Set labels and title
    control_labels=cell(1,size(waveform,1));
    for n=1:size(waveform,1)
        control_labels{n}=int2str(n);
    end

    % Make a translucent legend
    leg_obj=legend(control_labels,'Location','NorthEast');  
    set(leg_obj.BoxFace,'ColorType','truecoloralpha',...
                        'ColorData',uint8([200 200 200 64]')); 
    
    % Labels and the grid
    kxlabel('time, seconds'); ktitle('controls');
    kylabel('ens. average value, Hz'); kgrid; 

    % Report maximum nutation angle
    text(x_lower+0.05*(x_upper-x_lower),...
         y_lower+0.075*(y_upper-y_lower),...
         ['$\alpha_{\textnormal{max}}$ = ' ...
         num2str(max_nut,3) ' deg.'],...
         'Interpreter','latex');

    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Get phases and amplitudes
if ismember('phi_controls',spin_system.control.plotting)||...
   ismember('amp_controls',spin_system.control.plotting)
    
    % Preallocate phase and amplitude arrays
    amp_profile=zeros(size(waveform,1)/2,size(waveform,2));
    phi_profile=zeros(size(waveform,1)/2,size(waveform,2));
    
    % Fill the arrays
    for n=1:(size(waveform,1)/2)
        [amp_profile(n,:),phi_profile(n,:)]=cartesian2polar(waveform(2*n-1,:),waveform(2*n,:));
    end
    
end

% Plot phase controls
if ismember('phi_controls',spin_system.control.plotting)
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot);
    
    % Adapt to the integrator type
    switch spin_system.control.integrator
        
        % Piecewise-linear
        case 'trapezium'

            % Conventional linear plot
            p=plot(t_axis',wrapTo2Pi(phi_profile'));

        % Piecewise-constant
        case 'rectangle'
            
             % Stairs plot
            p=stairs(t_axis',wrapTo2Pi(phi_profile'));

        otherwise

            % Complain and bomb out
            error('unknown integrator type.');

    end

    % Plot with predictable colours
    for n=1:numel(p)
        p(n).Color=hsv2rgb([n/numel(p) 0.75 0.75]);
    end
       
    % Set labels and title
    control_labels=cell(1,size(waveform,1)/2);
    for k=1:(size(waveform,1)/2)
        control_labels{k}=['Ch ' int2str(2*k-1) ',' int2str(2*k)];
    end

    % Make a translucent legend
    leg_obj=legend(control_labels,'Location','NorthEast');  
    set(leg_obj.BoxFace,'ColorType','truecoloralpha',...
                        'ColorData',uint8([200 200 200 64]')); 
    
    % Labels and the grid
    kxlabel('time, seconds'); ktitle('control phases');
    kylabel('phase, radians'); kgrid;
    
    % Custom axis limits and ticks
    xlim([0 spin_system.control.pulse_dur]);
    ylim([0 2*pi]); yticks(0:pi/2:2*pi);
    yticklabels({'$0$','$0.5\pi$','$\pi$','$1.5\pi$','$2\pi$'});

    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Plot amplitude controls
if ismember('amp_controls',spin_system.control.plotting)
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot);
    
    % Adapt to the integrator type
    switch spin_system.control.integrator
        
        % Piecewise-linear
        case 'trapezium'

            % Conventional linear plot
            hold off; p=plot(t_axis',amp_profile'/(2*pi));

        % Piecewise-constant
        case 'rectangle'
            
            % Stairs plot
            hold off; p=stairs(t_axis',amp_profile'/(2*pi));

        otherwise

            % Complain and bomb out
            error('unknown integrator type.');

    end

    % Plot with predictable colours
    for n=1:numel(p)
        p(n).Color=hsv2rgb([n/numel(p) 0.75 0.75]);
    end
    
    % Determine the amplitude bound
    amp_bound=max([sqrt(2)*abs(upper_bound) sqrt(2)*abs(lower_bound)]);
    
    % Power bound
    max_amp=refline([0 amp_bound/(2*pi)]);
    set(max_amp,'Color','k','LineStyle','--');
    
    % Axis limits
    set(gca,'xlim',[0 spin_system.control.pulse_dur],...
            'ylim',[0 max([max(amp_profile(:)) 1.05*amp_bound])]/(2*pi));
    
    % Set labels and title
    control_labels=cell(1,size(waveform,1)/2);
    for k=1:(size(waveform,1)/2)
        control_labels{k}=['Ch ' int2str(2*k-1) ',' int2str(2*k)];
    end
    
    % Make a translucent legend
    leg_obj=legend(control_labels,'Location','NorthEast');  
    set(leg_obj.BoxFace,'ColorType','truecoloralpha',...
                        'ColorData',uint8([200 200 200 64]')); 
    
    % Labels and the grid
    kxlabel('time, seconds'); ktitle('control moduli');
    kylabel('ens. average value, Hz'); kgrid; 
    
    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Plot correlation orders
if ismember('correlation_order',spin_system.control.plotting)

    % Disallow Zeeman formalisms
    if ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
        error('trajectory analysis is not available in Zeeman formalisms, use sphten-liouv.');
    end
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot); hold off;
    
    % Loop over the ensemble
    for n=1:numel(traj_data)
        
        % Pull out a trajectory
        trajectory=traj_data{n}.forward;
        
        % Trace over the spatial degrees of freedom
        spn_dim=size(spin_system.bas.basis,1);
        spc_dim=size(trajectory,1)/spn_dim;
        trajectory=fpl2rho(trajectory,spc_dim);
        
        % Call trajectory analysis
        trajan(spin_system,trajectory,'correlation_order',t_axis); hold on;
        
    end
    
    % Update the axes
    kxlabel('time, seconds');
    xlim tight; ylim padded;
    
    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Plot coherence orders
if ismember('coherence_order',spin_system.control.plotting)

    % Disallow Zeeman formalisms
    if ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
        error('trajectory analysis is not available in Zeeman formalisms, use sphten-liouv.');
    end
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot); hold off;
    
    % Loop over the ensemble
    for n=1:numel(traj_data)
        
        % Pull out a trajectory
        trajectory=traj_data{n}.forward;
        
        % Trace over the spatial degrees of freedom
        spn_dim=size(spin_system.bas.basis,1);
        spc_dim=size(trajectory,1)/spn_dim;
        trajectory=fpl2rho(trajectory,spc_dim);
        
        % Call trajectory analysis
        trajan(spin_system,trajectory,'coherence_order',t_axis); hold on;
        
    end
    
    % Update the axes
    kxlabel('time, seconds');
    xlim tight; ylim padded;
    
    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Plot local probability densities
if ismember('local_each_spin',spin_system.control.plotting)

    % Disallow Zeeman formalisms
    if ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
        error('trajectory analysis is not available in Zeeman formalisms, use sphten-liouv.');
    end
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot); hold off;
    
    % Loop over the ensemble
    for n=1:numel(traj_data)
        
        % Pull out a trajectory
        trajectory=traj_data{n}.forward;
        
        % Trace over the spatial degrees of freedom
        spn_dim=size(spin_system.bas.basis,1);
        spc_dim=size(trajectory,1)/spn_dim;
        trajectory=fpl2rho(trajectory,spc_dim);
        
        % Call trajectory analysis
        trajan(spin_system,trajectory,'local_each_spin',t_axis); hold on;
        
    end
    
    % Update the axes
    kxlabel('time, seconds');
    xlim tight; ylim padded;
    
    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Plot total probability densities
if ismember('total_each_spin',spin_system.control.plotting)

    % Disallow Zeeman formalisms
    if ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
        error('trajectory analysis is not available in Zeeman formalisms, use sphten-liouv.');
    end
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot); hold off;
    
    % Loop over the ensemble
    for n=1:numel(traj_data)
        
        % Pull out a trajectory
        trajectory=traj_data{n}.forward;
        
        % Trace over the spatial degrees of freedom
        spn_dim=size(spin_system.bas.basis,1);
        spc_dim=size(trajectory,1)/spn_dim;
        trajectory=fpl2rho(trajectory,spc_dim);
        
        % Call trajectory analysis
        trajan(spin_system,trajectory,'total_each_spin',t_axis); hold on;
        
    end
    
    % Update the axes
    kxlabel('time, seconds');
    xlim tight; ylim padded;
    
    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Plot populations of the Zeeman energy levels
if ismember('level_populations',spin_system.control.plotting)

    % Disallow Zeeman formalisms
    if ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
        error('trajectory analysis is not available in Zeeman formalisms, use sphten-liouv.');
    end
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot); hold off;
    
    % Loop over the ensemble
    for n=1:numel(traj_data)
        
        % Pull out a trajectory
        trajectory=traj_data{n}.forward;
        
        % Trace over the spatial degrees of freedom
        spn_dim=size(spin_system.bas.basis,1);
        spc_dim=size(trajectory,1)/spn_dim;
        trajectory=fpl2rho(trajectory,spc_dim);
        
        % Call trajectory analysis
        trajan(spin_system,trajectory,'level_populations',t_axis); hold on;
        
    end
    
    % Update the axes
    kxlabel('time, seconds');
    xlim tight; ylim padded;
    
    % Increment plot counter
    current_plot=current_plot+1;
    
end

% Plot robustness
if ismember('robustness',spin_system.control.plotting)
    
    % Set the current plot
    subplot(n_plots_x,n_plots_y,current_plot);
    
    % Plot fidelity histogram
    histogram(fidelities(:),ceil(sqrt(numel(fidelities))));
    ktitle('ensemble'); kxlabel('fidelity'); 
    kylabel('number of systems'); kgrid;

    % Get histogram range + 10% margin
    edge_val=1.1*max(abs(fidelities),[],'all');
    xlim([-edge_val, +edge_val]); 
    
end

% Flush the graphics
drawnow();

end

% Consistency enforcement
function grumble(spin_system,waveform,traj_data,fidelities)
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be a real numerical array.');
end
if ~iscell(traj_data)
    error('traj_data must be a cell array of data structures.');
end
if (~isnumeric(fidelities))||(~isreal(fidelities))
    error('fidelities must be a real numerical array.');
end
if ismember('spectrogram',spin_system.control.plotting)
    if numel(unique(diff(spin_system.control.pulse_dt)))~=1
        error('spectrograms require uniform time intervals.');
    end
    if ~(mod(size(waveform,1),2)==0)
        error('spectrograms require an even number of control channels.');
    end
end
end

% If you can't explain it simply, you don't
% understand it well enough.
%
% Albert Einstein

