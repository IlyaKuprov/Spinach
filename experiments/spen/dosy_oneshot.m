% One-shot DOSY pulse sequence.
% 
%           fid=dosy_oneshot(spin_system,parameters,H,R,K,G,F)
%
% Parameters:
%
%   parameters.rho0           initial state 
%  
%   parameters.coil           detection state
%
%   parameters.spins          nuclei on which the sequence runs  
%
%   parameters.g_amp          gradient amplitude for diffusion 
%                             encoding (T/m)
%                         
%   parameters.g_dur          pulse width of the gradient for diffusion
%                             encoding (s)                    
%
%   parameters.kappa          unbalancing factor to unbalance the bipolar 
%                             gradients in the ratio (1+kappa):(1-kappa)
%
%   parameters.g_stab_del     gradient stabilization delay (s)
%
%   parameters.del            diffusion delay, seconds
% 
%   parameters.dims           size of the sample (m)
%
%   parameters.npts           number of discretization points in the grid
%
%   parameters.npoints        number of points in the acquired signal
% 
%   parameters.sweep          acquisition sweep width, Hz
%
%   H                         Fokker-Planck Hamiltonian
%
%   R                         Fokker-Planck relaxation superoperator
%
%   K                         Fokker-Planck kinetics superoperator
%
%   G                         Fokker-Planck gradient superoperators
%
%   F                         Fokker-Planck diffusion and flow 
%                             superoperator
%
% Outputs:
%
%   fid  - free induction decay
%           
% mariagrazia.concilio@sjtu.edu.cn
%
% <https://spindynamics.org/wiki/index.php?title=dosy_oneshot.m>

function fid=dosy_oneshot(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,F,G)

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Apply the first 90 degree pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Select coherence -1
rho=coherence(spin_system,rho,{{parameters.spins{1},-1}});

% Evolve under a positive gradient
rho=evolution(spin_system,L+(1+parameters.kappa)*parameters.g_amp*G{1},...
                          [],rho,parameters.g_dur/2,1,'final');

% Evolve under the first gradient stabilization delay
rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final');

% Apply the second 180 degree pulse
rho=step(spin_system,Ly,rho,pi);

% Select coherence +1
rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

% Evolve under a negative gradient
rho=evolution(spin_system,L-(1-parameters.kappa)*parameters.g_amp*G{1},[],...
                          rho,parameters.g_dur/2,1,'final');

% Evolve under the second gradient stabilization delay
rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final');

% Apply the third 90 degree pulse
rho=step(spin_system,Ly,rho,pi/2);

% Select coherence 0
rho=coherence(spin_system,rho,{{parameters.spins{1},0}});

% Evolve under a positive gradient 
rho=evolution(spin_system,L-(2*parameters.kappa)*parameters.g_amp*G{1},[],...
                          rho,parameters.g_dur/2,1,'final');

% Evolve under the second gradient stabilization delay
rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final');

% Run the diffusion time evolution
rho=evolution(spin_system,L,[],rho,parameters.del-4*(parameters.g_dur/2)...
                                                 -4*parameters.g_stab_del,1,'final');

% Evolve under a negative gradient
rho=evolution(spin_system,L-(2*parameters.kappa)*parameters.g_amp*G{1},[],...
                          rho,parameters.g_dur/2,1,'final');

% Evolve under the second gradient stabilization delay
rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final');

% Apply the second 90 degree pulse
rho=step(spin_system,Ly,rho,pi/2);

% Evolve under a positive gradient
rho=evolution(spin_system,L+(1+parameters.kappa)*parameters.g_amp*G{1},[],...
                          rho,parameters.g_dur/2,1,'final');

% Evolve under the third gradient stabilization delay
rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final');

% Apply the fifth 180 degree pulse
rho=step(spin_system,Ly,rho,pi);

% Evolve under a negative gradient
rho=evolution(spin_system,L-(1-parameters.kappa)*parameters.g_amp*G{1},[],...
                          rho,parameters.g_dur/2,1,'final');

% Evolve under the third gradient stabilization delay
rho=evolution(spin_system,L,[],rho,parameters.g_stab_del,1,'final');

% Run the F2 evolution
fid=evolution(spin_system,L,parameters.coil,rho,1/parameters.sweep,...
                            parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,F,G)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~isnumeric(F))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))||(~ismatrix(F))
    error('H, R, K and F must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))||(~all(size(K)==size(F)))
    error('H, R, K and F matrices must have the same dimension.');
end
if ~iscell(G)
    error('the G must be a 1x3 cell array.');
end
if ~isfield(parameters,'rho0')
    error('initial state should be specified in parameters.rho0.');
end
if ~isfield(parameters,'coil')
    error('detection state should be specified in parameters.coil.');
end
if ~isfield(parameters,'g_amp')
    error('the gradient amplitude should be specified in parameters.g_amp.');
elseif (~isnumeric(parameters.g_amp))||(~isreal(parameters.g_amp))||(numel(parameters.g_amp)~=1)
    error('parameters.g_amp array should have exactly one real element.');
end
if ~isfield(parameters,'g_dur')
    error('the gradient pulse width should be specified in parameters.g_dur.');
elseif (~isnumeric(parameters.g_dur))||(~isreal(parameters.g_dur))||...
       (numel(parameters.g_dur)~=1)||(parameters.g_dur<=0)
    error('parameters.g_dur should be a positive real scalar.');
end
if ~isfield(parameters,'kappa')
    error('the parameter kappa should be specified in parameters.kappa.');
elseif (~isnumeric(parameters.kappa))||(~isreal(parameters.kappa))||(numel(parameters.kappa)~=1)
    error('parameters.kappa array should have exactly one real element.');
end
if ~isfield(parameters,'g_stab_del')
    error('the gradient stabilization delay should be specified in parameters.g_stab_del.');
elseif (~isnumeric(parameters.g_stab_del))||(~isreal(parameters.g_stab_del))||...
       (numel(parameters.g_stab_del)~=1)||(parameters.g_stab_del<0)
    error('parameters.g_stab_del should be a non-negative real scalar.');
end
if ~isfield(parameters,'del')
    error('the diffusion delay should be specified in parameters.del.');
elseif (~isnumeric(parameters.del))||(~isreal(parameters.del))||...
       (numel(parameters.del)~=1)||(parameters.del<=0)
    error('parameters.del should be a positive real scalar.');
end
if parameters.del<(2*parameters.g_dur+4*parameters.g_stab_del)
    error('parameters.del should exceed 2*parameters.g_dur+4*parameters.g_stab_del.');
end
if ~isfield(parameters,'dims')
    error('sample dimension should be specified in parameters.dims.');
elseif numel(parameters.dims)~=1
    error('parameters.dims array should have exactly one element.');
end
if ~isfield(parameters,'npts')
    error('number of spin packets should be specified in parameters.npts.');
elseif numel(parameters.npts)~=1
    error('parameters.npts array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'sweep')
    error('the spectral range should be specified in parameters.sweep.');
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       (numel(parameters.sweep)~=1)||(parameters.sweep<=0)
    error('parameters.sweep should be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       (numel(parameters.npoints)~=1)||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
end

% Dr Who:  "I noticed you"
% Billie:  "But why?"
% Dr Who:  "Well, most people, when they don't 
%           understand something, they frown.
%           You - smile."


