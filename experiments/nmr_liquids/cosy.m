% Phase-sensitive COSY pulse sequence. Syntax:
%
%               fid=cosy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep        sweep width in Hz
%
%    parameters.npoints      number of points for both dimensions
%
%    parameters.spins        nuclei on which the sequence runs,
%                            specified as {'1H'}, {'13C'}, etc.
%
%    parameters.angle        second pulse angle in radians, usu-
%                            ally pi/2, but also allows COSY45, 
%                            COSY60, etc.
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid  - two-dimensional free induction decay
%
% Note: if the second pulse angle differs from 90 degrees, magni-
%       tude mode plotting is advised.
%
% gareth.charnock@oerc.ox.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cosy.m>

function fid=cosy(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Initial state
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');

% Detection state
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Get the pulse operator
Lp=operator(spin_system,'L+',parameters.spins{1});

% Apply the first pulse
rho=step(spin_system,(Lp+Lp')/2,rho,pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep,parameters.npoints(1)-1,'trajectory');

% Select "+1" coherence
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}});

% Apply the second pulse
rho_stack=step(spin_system,(Lp+Lp')/2,rho_stack,parameters.angle);
 
% Run the F2 evolution
fid=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep,parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('parameters.sweep array should have exactly one element.');
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
if ~isfield(parameters,'angle')
    error('first pulse angle should be specified in parameters.angle variable.');
elseif numel(parameters.angle)~=1
    error('parameters.angle array should have exactly one element.');
end
end

% Bring me into the company of those who seek truth, and deliver me from
% those who have found it.
%
% Arthur C. Clarke

