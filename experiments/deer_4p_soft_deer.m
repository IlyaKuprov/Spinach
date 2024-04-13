% Four-pulse DEER/PELDOR pulse sequence. The sequence uses soft 
% pulses computed with the Fokker-Planck formalism. Syntax:
%
%  echo_stack=deer_4p_soft_deer(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.pulse_frq  - frequencies for the four 
%                            pulses, Hz
%
%    parameters.pulse_pwr  - power levels for the four
%                            pulses, rad/s
%
%    parameters.pulse_dur  - durations for the four
%                            pulses, seconds
%
%    parameters.pulse_phi  - initial phases for the four 
%                            pulses, radians
%
%    parameters.pulse_rnk  - Fokker-Planck ranks for the
%                            four pulses
%
%    parameters.p1_p2_gap  - time between the end of the 
%                            first and the start of the
%                             second pulse, seconds
%
%    parameters.p2_p4_gap  - time between the end of the 
%                            second the start of the third
%                            pulse, seconds
%
%    parameters.p3_nsteps  - number of third pulse posi-
%                            tions in the interval between
%                            the first echo and the fourth
%                            pulse
%
%    parameters.echo_time  - time to sample around the ex-
%                            pected second echo position
%
%    parameters.echo_npts  - number of points in the second
%                            echo discretization
%
%    parameters.rho0       - initial state
%
%    parameters.coil       - detection state
%
%    parameters.method     - soft puse propagation method,
%                            'expv' for Krylov propagation,
%                            'expm' for exponential propa-
%                            gation, 'evolution' for Spin-
%                            ach evolution function
%
%   H  - Hamiltonian matrix, received from context function
%
%   R  - relaxation superoperator, received from context function
%
%   K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    echo_stack  - DEER echo stack, a matrix with p3_nsteps echoes
%                  with echo_npts points each
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
% Note: the time in the DEER trace refers to the third pulse inser-
%       tion point, after end of the second pulse.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=deer_4p_soft_deer.m>

function echo_stack=deer_4p_soft_deer(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Electron pulse operators
Ep=operator(spin_system,'L+',parameters.spins{1});
Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i;

% Frequency offsets
parameters.pulse_frq=-spin_system.inter.magnet*spin('E')/(2*pi)-...
                      parameters.pulse_frq-parameters.offset;

% First pulse
rho=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0,parameters.pulse_frq(1),parameters.pulse_pwr(1),...
                                                        parameters.pulse_dur(1),parameters.pulse_phi(1),...
                                                        parameters.pulse_rnk(1),parameters.method);
% Evolution
rho=evolution(spin_system,L,[],rho,parameters.p1_p2_gap,1,'final');

% Second pulse
rho=shaped_pulse_af(spin_system,L,Ex,Ey,rho,parameters.pulse_frq(2),parameters.pulse_pwr(2),...
                                            parameters.pulse_dur(2),parameters.pulse_phi(2),...
                                            parameters.pulse_rnk(2),parameters.method);
% Evolution
rho=evolution(spin_system,L,[],rho,parameters.p1_p2_gap+parameters.pulse_dur(1)/2,1,'final');

% Evolution
stepsize=(parameters.p2_p4_gap-parameters.p1_p2_gap)/parameters.p3_nsteps;
rho=evolution(spin_system,L,[],rho,stepsize,parameters.p3_nsteps,'trajectory');

% Third pulse
rho=shaped_pulse_af(spin_system,L,Ex,Ey,rho,parameters.pulse_frq(3),parameters.pulse_pwr(3),...
                                            parameters.pulse_dur(3),parameters.pulse_phi(3),...
                                            parameters.pulse_rnk(3),parameters.method);

% Evolution
rho(:,end:-1:1)=evolution(spin_system,L,[],rho(:,end:-1:1),stepsize,parameters.p3_nsteps,'refocus');

% Fourth pulse
rho=shaped_pulse_af(spin_system,L,Ex,Ey,rho,parameters.pulse_frq(4),parameters.pulse_pwr(4),...
                                            parameters.pulse_dur(4),parameters.pulse_phi(4),...
                                            parameters.pulse_rnk(4),parameters.method);
% Evolve to the start of the echo
echo_location=(parameters.p2_p4_gap-parameters.p1_p2_gap)+...
               parameters.pulse_dur(3)-parameters.echo_time/2;
rho=evolution(spin_system,L,[],rho,echo_location,1,'final');
                                                  
% Sample the echo
stepsize=parameters.echo_time/parameters.echo_npts;
echo_stack=evolution(spin_system,L,parameters.coil,rho,...
                     stepsize,parameters.echo_npts,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available in Liouville space.');
end
if ~isfield(parameters,'pulse_frq')
    error('pulse frequencies must be specified in parameters.pulse_frq field.');
end
if (~isnumeric(parameters.pulse_frq))||(~isreal(parameters.pulse_frq))||...
   (numel(parameters.pulse_frq)~=4)
    error('parameters.pulse_frq must have four real elements.');
end
if ~isfield(parameters,'pulse_pwr')
    error('pulse powers must be specified in parameters.pulse_pwr field.');
end
if (~isnumeric(parameters.pulse_pwr))||(~isreal(parameters.pulse_pwr))||...
   (numel(parameters.pulse_pwr)~=4)||any(parameters.pulse_pwr<=0)
    error('parameters.pulse_pwr must have four positive real elements.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse durations must be specified in parameters.pulse_dur field.');
end
if (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
   (numel(parameters.pulse_dur)~=4)||any(parameters.pulse_dur<=0)
    error('parameters.pulse_dur must have four positive real elements.');
end
if ~isfield(parameters,'pulse_phi')
    error('pulse phases must be specified in parameters.pulse_phi field.');
end
if (~isnumeric(parameters.pulse_phi))||(~isreal(parameters.pulse_phi))||...
   (numel(parameters.pulse_phi)~=4)
    error('parameters.pulse_phi must have four real elements.');
end
if ~isfield(parameters,'pulse_rnk')
    error('pulse grid ranks must be specified in parameters.pulse_rnk field.');
end
if (~isnumeric(parameters.pulse_rnk))||(~isreal(parameters.pulse_rnk))||...
   (numel(parameters.pulse_rnk)~=4)||any(mod(parameters.pulse_rnk,1)~=0)
    error('parameters.pulse_rnk must have four integer real elements.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'method')
    error('shaped pulse simulation method must be specified in parameters.method field.');
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
if ~isfield(parameters,'p1_p2_gap')
    error('p1-p2 time gap must be specified in parameters.p1_p2_gap field.');
end
if (~isnumeric(parameters.p1_p2_gap))||(~isreal(parameters.p1_p2_gap))||...
   (~isscalar(parameters.p1_p2_gap))||(parameters.p1_p2_gap<=0)
    error('parameters.p1_p2_gap must be a positive real scalar.');
end
if ~isfield(parameters,'p2_p4_gap')
    error('p2-p4 time gap must be specified in parameters.p2_p4_gap field.');
end
if (~isnumeric(parameters.p2_p4_gap))||(~isreal(parameters.p2_p4_gap))||...
   (~isscalar(parameters.p2_p4_gap))||(parameters.p2_p4_gap<=0)
    error('parameters.p2_p4_gap must be a positive real scalar.');
end
if ~isfield(parameters,'p3_nsteps')
    error('number of points in the trace must be specified in parameters.p3_nsteps field.');
end
if (~isnumeric(parameters.p3_nsteps))||(~isreal(parameters.p3_nsteps))||...
   (~isscalar(parameters.p3_nsteps))||(parameters.p3_nsteps<1)||...
   (mod(parameters.p3_nsteps,1)~=0)
    error('parameters.p3_nsteps must be a positive real integer.');
end
end

% The fate decides the job man best performs:
% One crafts a poem and another one informs,
% And he that has two talents in his purse
% Can choose to write denouncements in verse.
% 
% Stanislaw Jerzy Lec, translated by IK

