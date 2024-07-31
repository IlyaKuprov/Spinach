% Standard pulse-acquire sequence with a hard pulse. The user must sup-
% ply the pulse operator, the pulse duration and the initial condition.
% Echo detection may optionally be used. Syntax:
%
%             fid=hp_acquire(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep              sweep width, Hz
%
%    parameters.npoints            number of points in the FID
%
%    parameters.rho0               initial state
%
%    parameters.coil               detection state
%
%    parameters.pulse_op           pulse operator
%
%    parameters.pulse_angle        pulse angle, radians
%
%    parameters.decouple           spins to decouple, e.g. {'15N','13C'}
%                                  (sphten-liouv formalism only)
%
%    parameters.echo_time          optional echo time for echo detection
%                                  (echo_time - pulse - echo_time - fid)
%
%    parameters.echo_oper          optional pulse operator for echo
%                                  detection
%
%    parameters.echo_angle         optional pulse angle for echo detection
%
%    H     - Hamiltonian matrix, received from context function
%
%    R     - relaxation superoperator, received from context function
%
%    K     - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid   - free induction decay as seen by the state specified
%            in parameters parameters.coil
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=hp_acquire.m>

function fid=hp_acquire(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Project pulse operator
parameters.pulse_op=kron(speye(parameters.spc_dim),parameters.pulse_op);

% Apply the pulse
rho=step(spin_system,parameters.pulse_op,parameters.rho0,parameters.pulse_angle);

% Run the echo stage if specified
if isfield(parameters,'echo_time')
    
    % Run the first evolution time in the echo train
    rho=evolution(spin_system,L,[],rho,parameters.echo_time,1,'final');
    
    % Project pulse operator
    parameters.echo_oper=kron(speye(parameters.spc_dim),parameters.echo_oper);
    
    % Run the echo pulse
    rho=step(spin_system,parameters.echo_oper,rho,parameters.echo_angle);
    
    % Run the second evolution time in the echo train
    rho=evolution(spin_system,L,[],rho,parameters.echo_time,1,'final');
    
end

% Apply the decoupling
[L,rho]=decouple(spin_system,L,rho,parameters.decouple);

% Run the evolution and watch the coil state
fid=evolution(spin_system,L,parameters.coil,rho,...
              1/parameters.sweep,parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(numel(parameters.sweep)~=1)||(~isreal(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep should be a positive real number.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'pulse_op')
    error('pulse operator must be specified in parameters.pulse_op variable.');
end
if ~isfield(parameters,'pulse_angle')
    error('pulse angle must be specified in parameters.pulse_angle variable.');
end
if (~isnumeric(parameters.pulse_angle))||(numel(parameters.pulse_angle)~=1)||(~isreal(parameters.pulse_angle))
    error('parameters.pulse_angle should be a real number.');
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
    if ~ismember(spin_system.bas.formalism,{'zeeman-liouv'})
        error('analytical decoupling is only available for sphten-liouv formalism.');
    end
end
end

% Freedom (n.): To ask nothing. To expect nothing. To depend on nothing.
%
% Ayn Rand

