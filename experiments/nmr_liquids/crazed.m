% CRAZED pulse sequence. Ideal analytical coherence-pathway version
% of the sequence described in:
%
%              https://doi.org/10.1126/science.8266096
%
% Syntax:
%
%              fid=crazed(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep           sweep width in Hz
%
%    parameters.npoints         number of points for both dimensions
%
%    parameters.spins           nuclei on which the sequence runs,
%                               specified as {'1H'}, {'13C'}, etc.
%
%    parameters.angle           second pulse angle
%
%    parameters.rho0            initial condition
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
% Note: the gradient selection is represented by explicit coherence
%       projection onto the +2 double-quantum branch followed by the
%       +1 observable single-quantum branch.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=crazed.m>

function fid=crazed(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Compute the evolution timestep
timestep=1/parameters.sweep;

% Detection state
coil=state(spin_system,'L+',parameters.spins{1});

% Get the pulse operator
Ly=operator(spin_system,'Ly',parameters.spins{1});

% Apply the first pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep,...
                    parameters.npoints(1)-1,'trajectory');

% Apply a double-quantum filter:
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},2}});

% Apply the second pulse
rho_stack=step(spin_system,Ly,rho_stack,parameters.angle);

% Apply a single-quantum filter
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},1}});

% Run the F2 evolution
fid=evolution(spin_system,L,coil,rho_stack,timestep,...
              parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('parameters.sweep array should have exactly one element.');
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       (~isfinite(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a positive real scalar.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
elseif (~iscell(parameters.spins))||(~ischar(parameters.spins{1}))
    error('parameters.spins must be a one-element cell array of character strings.');
elseif ~ismember(parameters.spins{1},spin_system.comp.isotopes)
    error('parameters.spins refers to an isotope that is not present in the system.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two positive integers.');
end
if ~isfield(parameters,'angle')
    error('second pulse angle should be specified in parameters.angle variable.');
elseif numel(parameters.angle)~=1
    error('parameters.angle array should have exactly one element.');
elseif (~isnumeric(parameters.angle))||(~isreal(parameters.angle))||...
       (~isfinite(parameters.angle))
    error('parameters.angle must be a finite real scalar.');
end
if ~isfield(parameters,'rho0')
    error('initial state should be specified in parameters.rho0 variable.');
elseif ~isnumeric(parameters.rho0)
    error('parameters.rho0 must be a numeric array.');
elseif size(parameters.rho0,1)~=size(H,1)
    error('parameters.rho0 dimension must match the Liouville space dimension.');
end
end

% Цензуры у нас нет, но у каждого есть 
% чувство самосохранения.
%
% Виктор Коклюшкин

