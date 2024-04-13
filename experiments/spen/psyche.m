% PSYCHE pure-shift NMR pulse sequence. Syntax:
%
%        fid=psyche_1d(spin_system,parameters,H,R,K,G,F)
%
% The following parameters are required:
%
% parameters.rho0           initial state
%
% parameters.coil           detection state
%
% parameters.spins          nuclei on which the sequence runs
%
% parameters.g_amp          gradient amplitude (T/m)
%
% parameters.dims           size of the sample (m)
%
% parameters.npts           number of discretization points in the grid
%
% parameters.sweep          spectral range (Hz)
%
% parameters.npoints        number of points in the sweep
%
% parameters.zerofill       number of points for the zero filling
%
% parameters.diff           diffusion constant (m^2/s)
%
% H                         Fokker-Planck Hamiltonian, received
%                           from the imaging context
%
% R                         Fokker-Planck relaxation superoperator,
%                           received from the imaging context
%
% K                         Fokker-Planck kinetics superoperator,
%                           received from the imaging context
%
% G                         Fokker-Planck gradient superoperators,
%                           received from the imaging context
%
% F                         Fokker-Planck diffusion and flow super-
%                           operator, received from the context
%
% mohammadali.foroozandeh@chem.ox.ac.uk
% m.g.concilio@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=psyche.m>

function fid=psyche(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,F,G)

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Get chirp pulse waveform
[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.duration,...
                       parameters.bandwidth,parameters.smfactor,...
                       parameters.chirptype);

% Compute chirp RF amplitudes
if strcmp(parameters.chirptype,'saltire')
    rfbeta=(parameters.beta/360)*sqrt(2*parameters.bandwidth/(parameters.duration));
else
    q_beta=-(2*log(cosd(parameters.beta)/2+1/2))/pi;
    rfbeta=sqrt(parameters.duration*parameters.bandwidth*q_beta/(2*pi))/(parameters.duration);
end

% Calibrate chirps
norm_factor=max(Cx); 
Cx=Cx/norm_factor; Cy=Cy/norm_factor;
Cx=2*pi*rfbeta*Cx; Cy=2*pi*rfbeta*Cy;

% Apply the first pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);

% Run the first half of the t1 evolution
rho_stack=evolution(spin_system,L,[],rho,parameters.timestep1/2,...
                    parameters.npoints(1)-1,'trajectory');

% Select "+1" coherence
rho_stack=coherence(spin_system,rho_stack,{{'1H',+1}});

% Evolve for the delta period
rho_stack=evolution(spin_system,L+parameters.g_amp*G{1},[],...
                    rho_stack,parameters.delta,1,'final');

% Apply the hard 180 pulse
rho_stack=step(spin_system,Lx,rho_stack,pi);

% Evolve for the delta period
rho_stack=evolution(spin_system,L+parameters.g_amp*G{1},[],...
                    rho_stack,parameters.delta,1,'final');

% Select "-1" coherence
rho_stack=coherence(spin_system,rho_stack,{{'1H',-1}});

% Apply the 1st pulse of the PSYCHE element
durations=ones(size(Cx))*parameters.duration/numel(Cx);
rho_stack=shaped_pulse_xy(spin_system,L+parameters.g_amp*G{1},...
                          {Lx,Ly},{Cx,+Cy},durations,rho_stack,'expv-pwc');

% Run the diffusion time evolution
rho_stack=evolution(spin_system,L+parameters.g_amp*G{1},[],...
                    rho_stack,parameters.del-2*parameters.duration,1,'final');

% Select "0" coherence
rho_stack=coherence(spin_system,rho_stack,{{'1H',0}});

% Apply the 2nd pulse of the PSYCHE element
rho_stack=shaped_pulse_xy(spin_system,L+parameters.g_amp*G{1},...
                          {Lx,Ly},{Cx,-Cy},durations,rho_stack,'expv-pwc');

% Select "+1" coherence
rho_stack=coherence(spin_system,rho_stack,{{'1H',+1}});

% Run the second half of the t1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep1/2,...
                    parameters.npoints(1)-1,'refocus');

% Run the F2 evolution
fid=evolution(spin_system,L,parameters.coil,rho_stack,...
              parameters.timestep2,parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,F,G)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~isnumeric(F))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))||(~ismatrix(F))
    error('H, R, K and F must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))||(~all(size(K)==size(F)))
    error('H, R, K and F matrices must have the same dimension.');
end
if ~iscell(G)
    error('the G must be a 1x3 cell array.');
end
if ~isfield(parameters,'g_amp')
    error('the gradient amplitude should be specified in parameters.g_amp variable.');
elseif numel(parameters.g_amp)~=1
    error('parameters.g_amp array should have exactly one element.');
end
if ~isfield(parameters,'dims')
    error('sample dimension should be specified in parameters.dims variable.');
elseif numel(parameters.dims)~=1
    error('parameters.dims array should have exactly one element.');
end
if ~isfield(parameters,'npts')
    error('number of spin packets should be specified in parameters.npts variable.');
elseif numel(parameters.npts)~=1
    error('parameters.npts array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'diff')
    error('the diffusion coefficient should specified in parameters.diff variable.');
elseif numel(parameters.diff)~=1
    error('parameters.diff array should have exactly one element.');
end
end

% The rate of progress is proportional to the risk encountered. The
% public at large may well be more risk averse than the individuals
% in our business, but to limit the progress in the name of elimina-
% ting risk is no virtue.
%
% Neil Armstrong

