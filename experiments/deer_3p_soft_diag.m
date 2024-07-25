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
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=deer_3p_soft_diag.m>

function deer_3p_soft_diag(spin_system,parameters)

% Check consistency
grumble(parameters);

% Run pulse diagnostics experiment
fids=powder(spin_system,@deer_3p_soft_hole,parameters,parameters.assumptions);

% Apodisation and Fourier transform
fid_a=apodization(fids(:,1),'crisp-1d'); spectrum_a=fftshift(fft(fid_a,parameters.zerofill));
fid_b=apodization(fids(:,2),'crisp-1d'); spectrum_b=fftshift(fft(fid_b,parameters.zerofill));
fid_c=apodization(fids(:,3),'crisp-1d'); spectrum_c=fftshift(fft(fid_c,parameters.zerofill));
fid_d=apodization(fids(:,4),'crisp-1d'); spectrum_d=fftshift(fft(fid_d,parameters.zerofill));

% Plotting
figure(); scale_figure([2.0 1.5]);
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
figure(); surf(deer_axis_2d,echo_axis_2d,real(echo_stack));
ktitle('echo stack, unphased'); kylabel('echo window, ns');
kxlabel('2nd pulse position, $\mu$s'); axis tight; kgrid;

% Extract and phase the echo modulation
[deer_echoes,deer_sigmas,deer_traces]=svd(echo_stack);
deer_echoes=deer_echoes*deer_sigmas/deer_traces(1);
deer_traces=deer_traces*deer_sigmas/deer_traces(1);

% Plot echo components
figure(); plot(echo_axis,real(deer_echoes(:,1:3)));
kylabel('echo, real channel'); kgrid;
kxlabel('echo window, ns'); axis tight;
klegend({'${\bf{u}}_1\cdot\sigma_1$',...
         '${\bf{u}}_2\cdot\sigma_2$',...
         '${\bf{u}}_3\cdot\sigma_3$'});
ktitle('principal components of the echo stack');
    
% Plot DEER components
figure(); plot(deer_axis,real(deer_traces(:,1:3)));
kylabel('echo, real channel'); axis tight;
kxlabel('2nd pulse insertion point, $\mu$s');  
klegend({'${\bf{v}}_1\cdot\sigma_1$',...
         '${\bf{v}}_2\cdot\sigma_2$',...
         '${\bf{v}}_3\cdot\sigma_3$'}); kgrid;
ktitle('principal components of the echo stack');

end

% Consistency enforcement
function grumble(parameters)
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

