% Phase-sensitive homonuclear ROESY pulse sequence, assuming ideal
% spin-lock. Syntax:
%
%             fid=roesy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.sweep        - a vector with sweep widths
%                               in F1 and F2 directions, Hz
%
%     parameters.npoints      - a vector with point count
%                               in F1 and F2 directions
%
%     parameters.spins        - nuclei on which the sequence 
%                               runs, e.g. {'1H'}
%
%     parameters.tmix         - mixing time, seconds
%
%     H  - Hamiltonian matrix, received from context function
%
%     R  - relaxation superoperator, received from context function
%
%     K  - kinetics superoperator, received from context function
%
% Outputs:
%
%     fid.cos, fid.sin  -  components of the free induction
%                          decay for hypercomplex processing
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=roesy.m>

function fid=roesy(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timestep
timestep=1./parameters.sweep;

% Detection state
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% First pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);

% F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep(1),parameters.npoints(1)-1,'trajectory');

% Analytical spin lock
rho_stack_cos=spinlock(spin_system,Lx,Ly,rho_stack,'Y');
rho_stack_sin=spinlock(spin_system,Lx,Ly,rho_stack,'X');

% Mixing time
rho_stack_cos=evolution(spin_system,1i*R+1i*K,[],rho_stack_cos,parameters.tmix,1,'final');
rho_stack_sin=evolution(spin_system,1i*R+1i*K,[],rho_stack_sin,parameters.tmix,1,'final');

% F2 evolution
fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable');
fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('parameters.sweep array should have exactly two elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
if ~isfield(parameters,'tmix')
    error('mixing time should be specified in parameters.tmix variable.');
elseif numel(parameters.tmix)~=1
    error('parameters.tmix array should have exactly one element.');
end
end

% If you have always believed that everyone should play by the
% same rules and be judged by the same standards, that would
% have gotten you labeled a radical 60 years ago, a liberal 30
% years ago and a racist today.
%
% Thomas Sowell

