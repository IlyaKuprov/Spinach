% Three-pulse DEER/PELDOR pulse sequence. The sequence uses soft 
% pulses computed with the Fokker-Planck formalism. Syntax:
%
%   echo_stack=deer_3p_soft_deer(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.pulse_frq  - frequencies for the three 
%                            pulses, Hz
%
%    parameters.pulse_pwr  - power levels for the three
%                            pulses, rad/s
%
%    parameters.pulse_dur  - durations for the three
%                            pulses, seconds
%
%    parameters.pulse_phi  - initial phases for the three 
%                            pulses, radians
%
%    parameters.pulse_rnk  - Fokker-Planck ranks for the
%                            three pulses
%
%    parameters.p1_p3_gap  - time between the first and the
%                            third pulses, seconds
%
%    parameters.p2_nsteps  - number of second pulse posi-
%                            tions in the interval between
%                            the first and the third pulse
%
%    parameters.echo_time  - time to sample around the ex-
%                            pected echo position
%
%    parameters.echo_npts  - number of points in the echo
%                            discretization
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
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    echo_stack  - DEER echo stack, a matrix with p2_nsteps echoes
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
% Note: the time in the DEER trace refers to the second pulse inser-
%       tion point, after end of the first pulse.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=deer_3p_soft_deer.m>

function echo_stack=deer_3p_soft_deer(spin_system,parameters,H,R,K)

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
stepsize=parameters.p1_p3_gap/parameters.p2_nsteps;
rho_stack=evolution(spin_system,L,[],rho,stepsize,parameters.p2_nsteps,'trajectory');

% Second pulse
rho_stack=shaped_pulse_af(spin_system,L,Ex,Ey,rho_stack,parameters.pulse_frq(2),parameters.pulse_pwr(2),...
                                                        parameters.pulse_dur(2),parameters.pulse_phi(2),...
                                                        parameters.pulse_rnk(2),parameters.method);
% Evolution
rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),stepsize,parameters.p2_nsteps,'refocus');

% Third pulse
rho_stack=shaped_pulse_af(spin_system,L,Ex,Ey,rho_stack,parameters.pulse_frq(3),parameters.pulse_pwr(3),...
                                                        parameters.pulse_dur(3),parameters.pulse_phi(3),...
                                                        parameters.pulse_rnk(3),parameters.method);
% Evolve to the edge of the echo window
echo_location=parameters.p1_p3_gap+parameters.pulse_dur(1)/2+...
              parameters.pulse_dur(2)-parameters.echo_time/2;
rho_stack=evolution(spin_system,L,[],rho_stack,echo_location,1,'final');
                                                  
% Detect the echo
stepsize=parameters.echo_time/parameters.echo_npts;
echo_stack=evolution(spin_system,L,parameters.coil,rho_stack,...
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
end

% The substance of this book, as it is expressed in the editor's preface, is
% that to measure "right" by the false philosophy of the Hebrew prophets and
% "weepful" Messiahs is madness. Right is not the offspring of doctrine, but
% of power. All laws, commandments, or doctrines as to not doing to another
% what you do not wish done to you, have no inherent authority whatever, but
% receive it only from the club, the gallows and the sword. A man truly free
% is under no obligation to obey any injunction, human or divine. [...] Men
% should not be bound by moral rules invented by their foes.
%
% Leo Tolstoy, about Ragnar Redbeard's "Might is Right"

