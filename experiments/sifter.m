% SIFTER pulse sequence. Syntax:
%
%             fid=sifter(spin_system,parameters,H,R,K)
%
% where H is the Hamiltonian matrix, R is the relaxation matrix
% and K is the chemical kinetics matrix.
%
% Parameters:
%
%    parameters.npoints      - number of points in time evolution
%
%    parameters.timestep     - simulation time step, seconds
%
%    parameters.rho0         - initial state
%
%    parameters.coil         - detection state
%
%    parameters.pulse_opx    - pulse operator in X phase
%
%    parameters.pulse_opy    - pulse operator in Y phase
%
% Outputs:
%
%    fid - a 2D free induction decay
%
% alice.bowen@chem.ox.ac.uk
% Asif Equbal <asifequbal3@gmail.com>
%
% <https://spindynamics.org/wiki/index.php?title=sifter.m>

function fid=sifter(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% First 90-degree pulse along X
rho=step(spin_system,parameters.pulse_opx,parameters.rho0,pi/2);

% First part of the t1 period
rho_stack=evolution(spin_system,L,[],rho,parameters.timestep,...
                      parameters.npoints/2-1,'trajectory');
                
% Second 180-degree pulse along +X
rho_stack=step(spin_system,parameters.pulse_opx,rho_stack,pi);

% Second part of the t1 period
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep,...
                    parameters.npoints/2-1,'refocus');

% Third 90-degree pulse along +Y
rho_stack=step(spin_system,parameters.pulse_opy,rho_stack,pi/2);

% First part of the t2 period
rho_stack=evolution(spin_system,L,[],rho_stack(:,end:-1:1),parameters.timestep,...
                     parameters.npoints/2-1,'refocus');

% Fourth 180-degree pulse along +X
rho_stack=step(spin_system,parameters.pulse_opx,rho_stack,pi);

% Second part of the t2 period, 2D detection mode
fid=evolution(spin_system,L,parameters.coil,rho_stack,parameters.timestep,...
                     parameters.npoints/2-1,'observable');

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'npoints')
    error('number of points in the first echo must be specified in parameters.npoints field.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (~isscalar(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a positive real integer.');
end
if ~isfield(parameters,'timestep')
    error('time step must be specified in parameters.timestep field.');
end
if (~isnumeric(parameters.timestep))||(~isreal(parameters.timestep))||...
   (~isscalar(parameters.timestep))||(parameters.timestep<=0)
    error('parameters.timestep must be a positive real scalar.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'pulse_opx')
    error('X pulse operator must be supplied in parameters.pulse_opx field.');
end
if ~isfield(parameters,'pulse_opy')
    error('Y pulse operator must be supplied in parameters.pulse_opy field.');
end
end

% It would not surprise me for an instant to discover that there is an
% evolutionary instinct which causes us to prefer a slightly nasty and
% authentic person to someone who implausibly pretends to be flawless. 
%
% Rory Sutherland

