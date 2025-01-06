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
%   parameters.dims           size of the sample (m)
%
%   parameters.npts           number of discretization points in the grid
% 
%   parameters.diff           diffusion constant (m^2/s)
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
if ~isfield(parameters,'g_amp')
    error('the gradient amplitude should be specified in parameters.g_amp.');
elseif numel(parameters.g_amp)~=1
    error('parameters.g_amp array should have exactly one element.');
end
if ~isfield(parameters,'g_dur')
    error('the gradient pulse width should be specified in parameters.g_dur.');
elseif numel(parameters.g_dur)~=1
    error('parameters.g_dur array should have exactly one element.');
end
if ~isfield(parameters,'kappa')
    error('the parameter kappa should be specified in parameters.kappa.');
elseif numel(parameters.g_amp)~=1
    error('parameters.kappa array should have exactly one element.');
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
elseif numel(parameters.sweep)~=1
    error('parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints.');
elseif numel(parameters.npoints)~=1
    error('parameters.npoints array should have exactly one element.');
end
if ~isfield(parameters,'diff')
    error('the diffusion coefficient should specified in parameters.diff.');
elseif numel(parameters.diff)~=1
    error('parameters.diff array should have exactly one element.');
end
end

% Dr Who:  "I noticed you"
% Billie:  "But why?"
% Dr Who:  "Well, most people, when they don't 
%           understand something, they frown.
%           You - smile."

