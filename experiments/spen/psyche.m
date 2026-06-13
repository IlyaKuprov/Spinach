% PSYCHE pure-shift NMR pulse sequence. Syntax:
%
%        fid=psyche_1d(spin_system,parameters,H,R,K,G,F)
%
% Parameters:
%
%      parameters.rho0           initial state
%
%      parameters.coil           detection state
%
%      parameters.spins          nuclei on which the sequence runs
%
%      parameters.g_amp          gradient amplitude (T/m)
%
%      parameters.dims           size of the sample (m)
%
%      parameters.npts           number of discretization points in the grid
%
%      parameters.sweep          spectral range (Hz)
%
%      parameters.npoints        number of points in the sweep
%
%      parameters.zerofill       number of points for the zero filling
%
%      parameters.diff           diffusion constant (m^2/s)
%
%      H                         Fokker-Planck Hamiltonian, received
%                                from the imaging context
%
%      R                         Fokker-Planck relaxation superoperator,
%                                received from the imaging context
%
%      K                         Fokker-Planck kinetics superoperator,
%                                received from the imaging context
%
%      G                         Fokker-Planck gradient superoperators,
%                                received from the imaging context
%
%      F                         Fokker-Planck diffusion and flow super-
%                                operator, received from the context
%
% Outputs:
%
%      fid - a PSYCHE free induction decay as a 2D array
%
% mohammadali.foroozandeh@chem.ox.ac.uk
% mariagrazia.concilio@sjtu.edu.cn
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
if numel(G)~=3
    error('G must contain three gradient operators.');
end
if (~isnumeric(G{1}))||(~isnumeric(G{2}))||(~isnumeric(G{3}))
    error('gradient operators must be numeric matrices.');
end
if ~isfield(parameters,'g_amp')
    error('the gradient amplitude should be specified in parameters.g_amp variable.');
elseif (~isnumeric(parameters.g_amp))||(~isreal(parameters.g_amp))||...
       (numel(parameters.g_amp)~=1)||(~isfinite(parameters.g_amp))
    error('parameters.g_amp must be a finite real scalar.');
end
if ~isfield(parameters,'dims')
    error('sample dimension should be specified in parameters.dims variable.');
elseif (~isnumeric(parameters.dims))||(~isreal(parameters.dims))||...
       (numel(parameters.dims)~=1)||(~isfinite(parameters.dims))||(parameters.dims<=0)
    error('parameters.dims must be a positive real scalar.');
end
if ~isfield(parameters,'npts')
    error('number of spin packets should be specified in parameters.npts variable.');
elseif (~isnumeric(parameters.npts))||(~isreal(parameters.npts))||...
       (numel(parameters.npts)~=1)||(~isfinite(parameters.npts))||...
       (parameters.npts<1)||(mod(parameters.npts,1)~=0)
    error('parameters.npts must be a positive integer scalar.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
       (~ischar(parameters.spins{1}))
    error('parameters.spins must be a one-element cell array of character strings.');
end
if ~isfield(parameters,'diff')
    error('the diffusion coefficient should specified in parameters.diff variable.');
elseif (~isnumeric(parameters.diff))||(~isreal(parameters.diff))||...
       (numel(parameters.diff)~=1)||(~isfinite(parameters.diff))||(parameters.diff<0)
    error('parameters.diff must be a non-negative real scalar.');
end
if ~isfield(parameters,'rho0')
    error('initial state should be specified in parameters.rho0 variable.');
elseif (~isnumeric(parameters.rho0))||(~iscolumn(parameters.rho0))||...
       (size(parameters.rho0,1)~=size(H,1))
    error('parameters.rho0 must be a state vector of the same dimension as H.');
end
if ~isfield(parameters,'coil')
    error('detection state should be specified in parameters.coil variable.');
elseif (~isnumeric(parameters.coil))||(~iscolumn(parameters.coil))||...
       (size(parameters.coil,1)~=size(H,1))
    error('parameters.coil must be a state vector of the same dimension as H.');
end
if ~isfield(parameters,'npoints')
    error('point counts should be specified in parameters.npoints variable.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       (numel(parameters.npoints)~=2)||any(~isfinite(parameters.npoints))||...
       any(parameters.npoints<2)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two integer elements greater than 1.');
end
if ~isfield(parameters,'pulsenpoints')
    error('chirp pulse point count should be specified in parameters.pulsenpoints variable.');
elseif (~isnumeric(parameters.pulsenpoints))||(~isreal(parameters.pulsenpoints))||...
       (numel(parameters.pulsenpoints)~=1)||(~isfinite(parameters.pulsenpoints))||...
       (parameters.pulsenpoints<1)||(mod(parameters.pulsenpoints,1)~=0)
    error('parameters.pulsenpoints must be a positive integer scalar.');
end
if ~isfield(parameters,'duration')
    error('pulse duration should be specified in parameters.duration variable.');
elseif (~isnumeric(parameters.duration))||(~isreal(parameters.duration))||...
       (numel(parameters.duration)~=1)||(~isfinite(parameters.duration))||(parameters.duration<=0)
    error('parameters.duration must be a positive real scalar.');
end
if ~isfield(parameters,'bandwidth')
    error('pulse bandwidth should be specified in parameters.bandwidth variable.');
elseif (~isnumeric(parameters.bandwidth))||(~isreal(parameters.bandwidth))||...
       (numel(parameters.bandwidth)~=1)||(~isfinite(parameters.bandwidth))||(parameters.bandwidth<=0)
    error('parameters.bandwidth must be a positive real scalar.');
end
if ~isfield(parameters,'smfactor')
    error('chirp smoothing factor should be specified in parameters.smfactor variable.');
elseif (~isnumeric(parameters.smfactor))||(~isreal(parameters.smfactor))||...
       (numel(parameters.smfactor)~=1)||(~isfinite(parameters.smfactor))
    error('parameters.smfactor must be a finite real scalar.');
end
if ~isfield(parameters,'chirptype')
    error('chirp type should be specified in parameters.chirptype variable.');
elseif (~ischar(parameters.chirptype))||(~ismember(parameters.chirptype,{'wurst','wurst-adaptive',...
                                                                        'smoothed','smoothed-adaptive',...
                                                                        'saltire','saltire-adaptive'}))
    error('parameters.chirptype must be a supported chirp pulse type.');
end
if ~isfield(parameters,'beta')
    error('flip angle should be specified in parameters.beta variable.');
elseif (~isnumeric(parameters.beta))||(~isreal(parameters.beta))||...
       (numel(parameters.beta)~=1)||(~isfinite(parameters.beta))
    error('parameters.beta must be a finite real scalar.');
end
if ~isfield(parameters,'timestep1')
    error('first evolution timestep should be specified in parameters.timestep1 variable.');
elseif (~isnumeric(parameters.timestep1))||(~isreal(parameters.timestep1))||...
       (numel(parameters.timestep1)~=1)||(~isfinite(parameters.timestep1))||(parameters.timestep1<=0)
    error('parameters.timestep1 must be a positive real scalar.');
end
if ~isfield(parameters,'timestep2')
    error('second evolution timestep should be specified in parameters.timestep2 variable.');
elseif (~isnumeric(parameters.timestep2))||(~isreal(parameters.timestep2))||...
       (numel(parameters.timestep2)~=1)||(~isfinite(parameters.timestep2))||(parameters.timestep2<=0)
    error('parameters.timestep2 must be a positive real scalar.');
end
if ~isfield(parameters,'delta')
    error('diffusion delay should be specified in parameters.delta variable.');
elseif (~isnumeric(parameters.delta))||(~isreal(parameters.delta))||...
       (numel(parameters.delta)~=1)||(~isfinite(parameters.delta))||(parameters.delta<0)
    error('parameters.delta must be a non-negative real scalar.');
end
end

% The rate of progress is proportional to the risk encountered. The
% public at large may well be more risk averse than the individuals
% in our business, but to limit the progress in the name of elimina-
% ting risk is no virtue.
%
% Neil Armstrong

