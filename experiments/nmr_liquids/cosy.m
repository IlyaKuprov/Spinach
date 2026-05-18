% Phase-sensitive COSY pulse sequence from:
%
%           https://doi.org/10.1063/1.432450
%           https://doi.org/10.1016/0022-2364(82)90279-7
%
% Syntax:
%
%            fid=cosy(spin_system,parameters,H,R,K)
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
% Note: the implementation analytically retains the +1 t1 coher-
%       ence order, corresponding to a single phase-sensitive
%       pathway.
%
% gareth.charnock@oerc.ox.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cosy.m>

function fid=cosy(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Initial state up to an overall multiplier
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');

% Detection state up to an overall multiplier
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Get the pulse operator
Lx=operator(spin_system,'Lx',parameters.spins{1});

% Apply the first pulse
rho=step(spin_system,Lx,rho,pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep,...
                    parameters.npoints(1)-1,'trajectory');

% Select "+1" coherence
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}});

% Apply the second pulse
rho_stack=step(spin_system,Lx,rho_stack,parameters.angle);
 
% Run the F2 evolution
fid=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep,...
              parameters.npoints(2)-1,'observable');

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
end

% Bring me into the company of those who seek truth, 
% and deliver me from those who have found it.
%
% Arthur C. Clarke

