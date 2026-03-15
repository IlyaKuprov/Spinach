% Complete set of simulations related to three-pulse DEER. Runs pulse
% diagnostics, which is followed by echo diagnostics, which is follow-
% ed by DEER simulation. Syntax:
%
%                deer_3p_soft_diag(spin_system,parameters)
%
% Parameters:
%
%      parameters.pulse_frq  - frequencies for the three 
%                              pulses, Hz
%
%      parameters.pulse_pwr  - power levels for the three
%                              pulses, Hz
%
%      parameters.pulse_dur  - durations for the three
%                              pulses, seconds
%
%      parameters.pulse_phi  - initial phases for the three 
%                              pulses, radians
%
%      parameters.pulse_rnk  - Fokker-Planck ranks for the
%                              three pulses
%
%      parameters.offset     - receiver offset for the time
%                              domain detection, Hz
%
%      parameters.sweep      - sweep width for time domain
%                              detection, Hz
%
%      parameters.npoints    - number of points in the free
%                              induction decay 
%
%      parameters.zerofill   - length of the zero-filled FFT
%
%      parameters.rho0       - initial state
%
%      parameters.coil       - detection state
%
%      parameters.p1_p3_gap  - time between the first and the
%                              third pulses, seconds
%
%      parameters.p2_nsteps  - number of second pulse posi-
%                              tions in the interval between
%                              the first and the third pulse
%
%      parameters.echo_time  - time to sample around the ex-
%                              pected echo position
%
%      parameters.echo_npts  - number of points in the echo
%                              discretization
%
%      parameters.method     - soft puse propagation method,
%                              'expv' for Krylov propagation,
%                              'expm' for exponential propa-
%                              gation, 'evolution' for Spin-
%                              ach evolution function
%
%      parameters.assumptions - Hamiltonian generation assump-
%                               tions, use 'deer' to keep two-
%                               electron flip-flop terms and 
%                               'deer-zz' to drop them
%
% Outputs:
%
%   Figure 1: pulse diagnostics
%   Figure 2: DEER echo stack
%   Figure 3: principal components of the stack, echo
%   Figure 4: principal components of the stack, DEER
%
% Note: for the method, start with 'expm', change to 'expv' if the
%       calculation runs out of memory, and use 'evolution' as the
%       last resort.
%
% Note: simulated echoes tend to be sharp and hard to catch becau-
%       se simulation does not have distributions in experimental
%       parameters. Fourier transforming the echo prior to integ-
%       ration is recommended.
%
% Note: the time in the DEER trace refers to the second pulse inser-
%       tion point, after end of first pulse.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=deer_3p_soft_diag.m>

function deer_3p_soft_diag(spin_system,parameters)

% Check consistency
grumble(parameters);

% Run pulse diagnostics experiment
fids=powder(spin_system,@deer_3p_soft_hole,parameters,parameters.assumptions);

% Apodisation and Fourier transform
fid_a=apodisation(spin_system,fids(:,1),{{'crisp'}}); spectrum_a=fftshift(fft(fid_a,parameters.zerofill));
fid_b=apodisation(spin_system,fids(:,2),{{'crisp'}}); spectrum_b=fftshift(fft(fid_b,parameters.zerofill));
fid_c=apodisation(spin_system,fids(:,3),{{'crisp'}}); spectrum_c=fftshift(fft(fid_c,parameters.zerofill));
fid_d=apodisation(spin_system,fids(:,4),{{'crisp'}}); spectrum_d=fftshift(fft(fid_d,parameters.zerofill));

% Plotting
kfigure(); scale_figure([2.0 1.5]);
subplot(2,2,1); plot_1d(spin_system,real(spectrum_a),parameters,'r-'); ktitle('frequency sweep epr');
subplot(2,2,2); plot_1d(spin_system,real(spectrum_b),parameters,'r-'); ktitle('first pulse');
subplot(2,2,3); plot_1d(spin_system,real(spectrum_c),parameters,'r-'); ktitle('second pulse');
subplot(2,2,4); plot_1d(spin_system,real(spectrum_d),parameters,'r-'); ktitle('third pulse');
drawnow();

% Get DEER echo stack
echo_stack=powder(spin_system,@deer_3p_soft_deer,parameters,parameters.assumptions);

% Time axis for the echo window
echo_axis=1e9*linspace(-parameters.echo_time/2,...
                        parameters.echo_time/2,parameters.echo_npts+1);

% Time axis for the DEER trace
deer_axis=1e6*linspace(0,parameters.p1_p3_gap,parameters.p2_nsteps+1);

% Axis set for the echo stack
[deer_axis_2d,echo_axis_2d]=meshgrid(deer_axis,echo_axis);

% Plot the echo stack
kfigure(); surf(deer_axis_2d,echo_axis_2d,real(echo_stack));
ktitle('echo stack, unphased'); kylabel('echo window, ns');
kxlabel('2nd pulse position, $\mu$s'); axis tight; kgrid;

% Extract and phase the echo modulation
[deer_echoes,deer_sigmas,deer_traces]=svd(echo_stack);
deer_echoes=deer_echoes*deer_sigmas/deer_traces(1);
deer_traces=deer_traces*deer_sigmas/deer_traces(1);

% Plot echo components
kfigure(); plot(echo_axis,real(deer_echoes(:,1:3)));
kylabel('echo, real channel'); kgrid;
kxlabel('echo window, ns'); axis tight;
klegend({'${\bf{u}}_1\cdot\sigma_1$',...
         '${\bf{u}}_2\cdot\sigma_2$',...
         '${\bf{u}}_3\cdot\sigma_3$'});
ktitle('principal components of the echo stack');
    
% Plot DEER components
kfigure(); plot(deer_axis,real(deer_traces(:,1:3)));
kylabel('echo, real channel'); axis tight;
kxlabel('2nd pulse insertion point, $\mu$s');  
klegend({'${\bf{v}}_1\cdot\sigma_1$',...
         '${\bf{v}}_2\cdot\sigma_2$',...
         '${\bf{v}}_3\cdot\sigma_3$'}); kgrid;
ktitle('principal components of the echo stack');

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'pulse_frq')
    error('pulse frequencies must be specified in parameters.pulse_frq field.');
end
if (~isnumeric(parameters.pulse_frq))||(~isreal(parameters.pulse_frq))||...
   (numel(parameters.pulse_frq)~=3)
    error('parameters.pulse_frq must have three real elements.');
end
if ~isfield(parameters,'pulse_pwr')
    error('pulse powers must be specified in parameters.pulse_pwr field.');
end
if (~isnumeric(parameters.pulse_pwr))||(~isreal(parameters.pulse_pwr))||...
   (numel(parameters.pulse_pwr)~=3)||any(parameters.pulse_pwr<=0)
    error('parameters.pulse_pwr must have three positive real elements.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse durations must be specified in parameters.pulse_dur field.');
end
if (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
   (numel(parameters.pulse_dur)~=3)||any(parameters.pulse_dur<=0)
    error('parameters.pulse_dur must have three positive real elements.');
end
if ~isfield(parameters,'pulse_phi')
    error('pulse phases must be specified in parameters.pulse_phi field.');
end
if (~isnumeric(parameters.pulse_phi))||(~isreal(parameters.pulse_phi))||...
   (numel(parameters.pulse_phi)~=3)
    error('parameters.pulse_phi must have three real elements.');
end
if ~isfield(parameters,'pulse_rnk')
    error('pulse grid ranks must be specified in parameters.pulse_rnk field.');
end
if (~isnumeric(parameters.pulse_rnk))||(~isreal(parameters.pulse_rnk))||...
   (numel(parameters.pulse_rnk)~=3)||any(mod(parameters.pulse_rnk,1)~=0)
    error('parameters.pulse_rnk must have three integer real elements.');
end
if ~isfield(parameters,'offset')
    error('receiver offset must be specified in parameters.offset field.');
end
if (~isnumeric(parameters.offset))||(~isreal(parameters.offset))||...
   (~isscalar(parameters.offset))
    error('parameters.offset must be a real scalar.');
end
if ~isfield(parameters,'sweep')
    error('width of the detection window must be specified in parameters.sweep field.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isscalar(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of points in the FID must be specified in parameters.npoints field.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (~isscalar(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a positive real integer.');
end
if ~isfield(parameters,'zerofill')
    error('FFT zero fill length must be specified in parameters.zerofill field.');
end
if (~isnumeric(parameters.zerofill))||(~isreal(parameters.zerofill))||...
   (~isscalar(parameters.zerofill))||(parameters.zerofill<parameters.npoints)||...
   (mod(parameters.zerofill,1)~=0)
    error('parameters.zerofill must be an integer not smaller than parameters.npoints.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'p1_p3_gap')
    error('p1-p3 time gap must be specified in parameters.p1_p3_gap field.');
end
if (~isnumeric(parameters.p1_p3_gap))||(~isreal(parameters.p1_p3_gap))||...
   (~isscalar(parameters.p1_p3_gap))||(parameters.p1_p3_gap<=0)
    error('parameters.p1_p3_gap must be a positive real scalar.');
end
if ~isfield(parameters,'p2_nsteps')
    error('number of points in the trace must be specified in parameters.p2_nsteps field.');
end
if (~isnumeric(parameters.p2_nsteps))||(~isreal(parameters.p2_nsteps))||...
   (~isscalar(parameters.p2_nsteps))||(parameters.p2_nsteps<1)||...
   (mod(parameters.p2_nsteps,1)~=0)
    error('parameters.p2_nsteps must be a positive real integer.');
end
if ~isfield(parameters,'echo_time')
    error('width of the echo window must be specified in parameters.echo_time field.');
end
if (~isnumeric(parameters.echo_time))||(~isreal(parameters.echo_time))||...
   (~isscalar(parameters.echo_time))||(parameters.echo_time<=0)
    error('parameters.echo_time must be a positive real scalar.');
end
if ~isfield(parameters,'echo_npts')
    error('number of points in the echo must be specified in parameters.echo_npts field.');
end
if (~isnumeric(parameters.echo_npts))||(~isreal(parameters.echo_npts))||...
   (~isscalar(parameters.echo_npts))||(parameters.echo_npts<1)||...
   (mod(parameters.echo_npts,1)~=0)
    error('parameters.echo_npts must be a positive real integer.');
end
if ~isfield(parameters,'method')
    error('shaped pulse simulation method must be specified in parameters.method field.');
end
if ~ischar(parameters.method)
    error('parameters.method must be a character string.');
end
if ~isfield(parameters,'assumptions')
    error('assumptions must be specified in parameters.assumptions field.');
end
if ~ischar(parameters.assumptions)
    error('parameters.assumptions must be a character string.');
end
end

% The "name-blind" application procedures (wherein the applicant's name
% is removed from the CV) were hastily discontinued in Australian Public
% Service after it emerged that gender and ethnicity biases in recruit-
% ment had increased, rather than decreased, as a result.

