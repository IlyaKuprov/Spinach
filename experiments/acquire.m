% Simple forward time evolution with signal acquisition. Syntax:
%
%            fid=acquire(spin_system,parameters,H,R,K)
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
%    parameters.decouple           spins to decouple, e.g. {'15N','13C'}
%
%    parameters.homodec_oper       operator to add to the Liouvillian at
%                                  the detection stage
%
%    parameters.homodec_pwr        power coefficient for the operator, Hz
%
%    parameters.dead_time          the system will be evolved for this
%                                  time (seconds) before the signal 
%                                  acquisition begins
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
% <https://spindynamics.org/wiki/index.php?title=acquire.m>

function fid=acquire(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Apply the decoupling
[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple);

% Apply the homodecoupling
if isfield(parameters,'homodec_oper')
    
    % Project into FP space
    parameters.homodec_oper=kron(speye(parameters.spc_dim),parameters.homodec_oper);
    
    % Add to the Liouvillian
    L=L+2*pi*parameters.homodec_pwr*parameters.homodec_oper;
    
end

% Run the dead time
if isfield(parameters,'dead_time')
    parameters.rho0=step(spin_system,L,parameters.rho0,parameters.dead_time);
end

% Run the evolution and watch the coil state
fid=evolution(spin_system,L,parameters.coil,parameters.rho0,...
              1/parameters.sweep,parameters.npoints-1,'observable');

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
if (~isnumeric(parameters.npoints))||(~isscalar(parameters.npoints))||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
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

% According to Oxford Chemistry folklore, Steve Davies's car, usually
% parked preposterously right in front of PTCL entrance door, had suf-
% fered more than just a few bird droppings in its life. David Logan
% got tired one day of the thing blocking the footpath and made good
% use of the car's low clearance (which made it level with the porch)
% by walking right over it, to much cheering and jubilation from the
% accompanying students. The car has been parked well to the side ev-
% er since.

