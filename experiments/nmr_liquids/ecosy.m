% Phase-sensitive E.COSY pulse sequence. Syntax:
%
%               fid=ecosy(spin_system,parameters,H,R,K)
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
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid.cos, fid.sin - real and imaginary components of the
%                       States quadrature signal
%
% Notes: with B0=B1=0, the multiple-quantum filter weights are
%        B_p=p^2/4 for even p and B_p=(p^2-1)/4 for odd p.
%        This gives 1, 2, 4, 6, and 9 for p=2..6; coherence
%        orders above six are a documented implementation limit.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ecosy.m>

function fid=ecosy(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Compute evolution timestep
timestep=1/parameters.sweep;

% Initial (post-pulse) and detection states
rho0=state(spin_system,'Lx',parameters.spins{1});
coil=state(spin_system,'L+',parameters.spins{1});

% Get pulse operators
Lx=operator(spin_system,'Lx',parameters.spins{1});
Ly=operator(spin_system,'Ly',parameters.spins{1});

% Run F1 evolution
rho_stack=evolution(spin_system,L,[],rho0,timestep,...
                    parameters.npoints(1)-1,'trajectory');

% Apply the second pulse (States quadrature)
rho_stack_cos=step(spin_system,Lx,rho_stack,pi/2);
rho_stack_sin=step(spin_system,Ly,rho_stack,pi/2);

% Apply the multiple-quantum filter
rho_stack_cos=1*coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+2,-2]}})+...
              2*coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+3,-3]}})+...
              4*coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+4,-4]}})+...
              6*coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+5,-5]}})+...
              9*coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+6,-6]}});
rho_stack_sin=1*coherence(spin_system,rho_stack_sin,{{parameters.spins{1},[+2,-2]}})+...
              2*coherence(spin_system,rho_stack_sin,{{parameters.spins{1},[+3,-3]}})+...
              4*coherence(spin_system,rho_stack_sin,{{parameters.spins{1},[+4,-4]}})+...
              6*coherence(spin_system,rho_stack_sin,{{parameters.spins{1},[+5,-5]}})+...
              9*coherence(spin_system,rho_stack_sin,{{parameters.spins{1},[+6,-6]}});

% Apply the third pulse
rho_stack_cos=step(spin_system,Lx,rho_stack_cos,pi/2);
rho_stack_sin=step(spin_system,Ly,rho_stack_sin,pi/2);

% Run the F2 evolution
fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep,...
                  parameters.npoints(2)-1,'observable');
fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep,...
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
elseif (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
       (~ischar(parameters.spins{1}))
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
end

% One man's caching policy is another man's memory leak.

