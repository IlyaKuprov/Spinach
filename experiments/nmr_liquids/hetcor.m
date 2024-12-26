% Magnitude-mode HETCOR pulse sequence. Syntax:
%
%            fid=hetcor(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep         [F1 F2] sweep widths in the 
%                             two frequency directions, Hz
%
%    parameters.npoints       [F1 F2] numbers of points in 
%                             the two time directions in fid
%
%    parameters.spins         {F1 F2} nuclei (e.g. {'1H','13C'})
%
%    parameters.decouple      list of nuclei that detection 
%                             time decoupling should be applied
%                             to - cell array of strings, e.g.
%                             {'1H','15N'})
%
%    parameters.J             working scalar coupling, Hz
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid - two-dimensional free induction decay for magnitude-
%          mode processing
%
% Note: natural abundance experiments should make use of the iso-
%       tope dilution functionality. See dilute.m function.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hetcor.m>

function fid=hetcor(spin_system,parameters,H,R,K)

% Consistenchy check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Evolution timesteps
timestep=1./parameters.sweep;

% J-coupling evolution times
delta_2=abs(1/(2*parameters.J));
delta_3=abs(1/(3*parameters.J));

% Initial state
rho=state(spin_system,'Lz',parameters.spins{1});

% Detection state
coil=state(spin_system,'L+',parameters.spins{2});

% Pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx_F1=(Lp+Lp')/2; Ly_F1=(Lp-Lp')/2i;
Lp=operator(spin_system,'L+',parameters.spins{2});
Lx_F2=(Lp+Lp')/2;

% First pulse on F1 with coherence selection
rho=   step(spin_system,Lx_F1,rho,pi/2)+...
    1i*step(spin_system,Ly_F1,rho,pi/2);

% F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2,...
                    parameters.npoints(1)-1,'trajectory');
rho_stack=step(spin_system,Lx_F2,rho_stack,pi);
rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2,...
                    parameters.npoints(1)-1,'refocus');

% 1/2J part of delta evolution
rho_stack=evolution(spin_system,L,[],rho_stack,delta_2,1,'final');

% Two pi/2 pulses with coherence selection
rho_stack=+step(spin_system,Lx_F1+Lx_F2,rho_stack,+pi/2)...
          +step(spin_system,Lx_F1+Lx_F2,rho_stack,-pi/2);

% 1/3J part of delta evolution
rho_stack=evolution(spin_system,L,[],rho_stack,delta_3,1,'final');

% F1 decoupling
[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.decouple);

% F2 detection
fid=evolution(spin_system,L,coil,rho_stack,timestep(2),...
              parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
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
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly two elements.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
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
if ~isfield(parameters,'J')
    error('scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('parameters.J array should have exactly one element.');
end
end

% Studies have suggested that hypomania can heighten certain cognitive pro-
% cesses, increase original and idiosyncratic thought, and even enhance lin-
% guistic skills. Manic states can also reduce the need for sleep, foster
% intense and obsessive concentration, create unmitigated self-confidence,
% and eliminate concern for social norms - just what you need, perhaps, to
% push the envelope of artistic creativity.
%
% M.F. Bear, B.W. Connors, M.A. Paradiso,
% "Neuroscience: exploring the brain", Kluwer, 2001.

