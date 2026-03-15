% Amplitude-mode homonuclear TOCSY pulse sequence. Syntax:
%
%             fid=tocsy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%      parameters.sweep   -  sweep widths, Hz
%
%      parameters.npoints -  number of points for both 
%                            dimensions
%
%      parameters.spins   -  nuclei on which the sequence runs,
%                            specified as '1H', '13C', etc.
%
%      parameters.tmix    -  mixing time, seconds
%
%      parameters.lamp    -  spin-lock power, Hz
%
%      parameters.rho0    -  initial state
%
%      H - Hamiltonian matrix, received from context function
%
%      R - relaxation superoperator, received from context function
%
%      K - kinetics superoperator, received from context function
%
% Outputs:
%
%      fid.cos, fid.sin   -  sine and cosine components
%                            of the States quadrature
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=tocsy.m>

function fid=tocsy(spin_system,parameters,H,R,K)

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

% Mixing time under spin-lock
rho_stack_cos=evolution(spin_system,H+2*pi*parameters.lamp*Lx,[],rho_stack,parameters.tmix,1,'final');
rho_stack_sin=evolution(spin_system,H+2*pi*parameters.lamp*Ly,[],rho_stack,parameters.tmix,1,'final');

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
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       any(parameters.sweep<=0)
    error('parameters.sweep must contain two positive real numbers.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
elseif (~iscell(parameters.spins))||(~ischar(parameters.spins{1}))
    error('parameters.spins must be a single-element cell array of character strings.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two positive integers.');
end
if ~isfield(parameters,'tmix')
    error('mixing time should be specified in parameters.tmix variable.');
elseif numel(parameters.tmix)~=1
    error('parameters.tmix array should have exactly one element.');
elseif (~isnumeric(parameters.tmix))||(~isreal(parameters.tmix))||...
       (parameters.tmix<0)
    error('parameters.tmix must be a non-negative real scalar.');
end
if ~isfield(parameters,'lamp')
    error('spin-lock power should be specified in parameters.lamp variable.');
elseif numel(parameters.lamp)~=1
    error('parameters.lamp array should have exactly one element.');
elseif (~isnumeric(parameters.lamp))||(~isreal(parameters.lamp))
    error('parameters.lamp must be a real scalar.');
end
if ~isfield(parameters,'rho0')
    error('initial state should be specified in parameters.rho0 variable.');
elseif ~isnumeric(parameters.rho0)
    error('parameters.rho0 must be a numeric array.');
end
end

% There is a theory which states that if ever anyone discovers exactly what
% the Universe is for and why it is here, it will instantly disappear and
% be replaced by something even more bizarre and inexplicable. There is
% another theory mentioned, which states that this has already happened.
%
% Douglas Adams, "Hitchhiker's Guide to the Galaxy"

