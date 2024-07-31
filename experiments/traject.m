% Simple forward time evolution trajectory. Syntax:
%
%            traj=traject(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep       sweep width, Hz
%
%    parameters.npoints     number of points in the trajectory
%
%    parameters.rho0        initial state
%
%    parameters.decouple    spins to decouple, e.g. {'15N','13C'}
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    traj  - system trajectory, a bookshelf stack of state vectors
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=traject.m>

function traj=traject(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Apply the decoupling
[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple);

% Run the evolution and watch the coil state
traj=evolution(spin_system,L,[],parameters.rho0,1/parameters.sweep,...
                                parameters.npoints-1,'trajectory');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(numel(parameters.sweep)~=1)||...
   (~isreal(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep should be a positive real number.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'decouple')
    error('list of decoupled spins (or an empty cell array) should be supplied in parameters.decouple variable.');
end
if ~iscell(parameters.decouple)
    error('parameters.decouple must be a cell array of strings.');
end
if numel(parameters.decouple)>0
    if any(~cellfun(@ischar,parameters.decouple))
        error('elements of parameters.decouple cell array must be strings.');
    end 
    if any(~ismember(parameters.decouple,spin_system.comp.isotopes))
        error('parameters.decouple contains isotopes that are not present in the system.');
    end
    if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
        error('analytical decoupling is only available for sphten-liouv formalism.');
    end
end
end

% That we are berated by our enemies means that we 
% are doing everything right.
%
% Joseph Stalin

