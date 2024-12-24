% Magnitude mode PANSY-COSY pulse sequence from Figure 1 in the 
% paper by Kupce, Freeman, and John:
%
%               http://dx.doi.org/10.1021/ja0634876
%
% Syntax:
%
%            fid=pansy_cosy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins   - nuclei on which the sequence runs,
%                         e.g. {'1H','13C'}
%
%    parameters.sweep   - a vector with two sweep widths in Hz
%
%    parameters.npoints - a vector of integers specifying 
%                         point count in each dimension
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid.aa - magnitude mode COSY FID on F1,F1 nuclei
%
%    fid.ab - magnitude mode COSY FID on F1,F2 nuclei
%
% Note: decoupling with respect to either of the working nuclei 
%       is impossible in this pulse sequence.
%
% Andrew Porter
% i.ya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=pansy_cosy.m>

function fid=pansy_cosy(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get evolution time steps
timesteps=1./parameters.sweep;

% Initial state
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');

% Detection state
coil_h=state(spin_system,'L+',parameters.spins{1},'cheap');
coil_c=state(spin_system,'L+',parameters.spins{2},'cheap');

% Get pulse operators
Hp=operator(spin_system,'L+',parameters.spins{1});
Cp=operator(spin_system,'L+',parameters.spins{2});
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; Cy=(Cp-Cp')/2i;

% Apply the first pulse
rho_p=step(spin_system,Hx,rho,+pi/2);
rho_m=step(spin_system,Hx,rho,-pi/2);

% Select "+1" coherence
rho_p=coherence(spin_system,rho_p,{{parameters.spins{1},+1}});
rho_m=coherence(spin_system,rho_m,{{parameters.spins{1},+1}});

% Run F1 evolution
rho_stack_p=evolution(spin_system,L,[],rho_p,timesteps(1),...
                      parameters.npoints(1)-1,'trajectory');
rho_stack_m=evolution(spin_system,L,[],rho_m,timesteps(1),...
                      parameters.npoints(1)-1,'trajectory');

% Apply the second pulse pair
rho_stack_p=step(spin_system,Hy+Cy,rho_stack_p,pi/2);
rho_stack_m=step(spin_system,Hy+Cy,rho_stack_m,pi/2);

% Perform axial peak elimination
rho_stack=rho_stack_p-rho_stack_m;

% Run F2 evolution
fid.aa=evolution(spin_system,L,coil_h,rho_stack,timesteps(1),...
                 parameters.npoints(1)-1,'observable');
fid.ab=evolution(spin_system,L,coil_c,rho_stack,timesteps(2),...
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
elseif numel(parameters.sweep)~=2
    error('parameters array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly two element.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.points variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
end

% The problem with being working class is that people only pay you to do
% things which are actually useful. You don't find scaffolders randomly
% erecting scaffolding where it isn't needed. Middle-class people, on the
% other hand, can easily generate their own bullshit. In an IT-packed
% office environment, lots of moronic activities look indistinguishable
% from productive work.
%
% Rory Sutherland

