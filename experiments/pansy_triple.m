% Triple-channel PANSY pulse sequence from Figure 1 in the 
% paper by Kupce, Freeman, and John:
%
%            https://dx.doi.org/10.1021/ja0634876
%
% Syntax:
%
%       fid=pansy_triple(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins   - three nuclei on which the sequence 
%                         runs, e.g. {'1H','13C','15N'}
%
%    parameters.sweep   - a vector with three sweep widths
%                         in Hz
%
%    parameters.npoints - a vector of three integers specify-
%                         ing point count in each dimension
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid.aa - COSY FID on F1 nuclei
%
%    fid.ab - COSY FID on F1,F2 nuclei
%
%    fid.ac - COSY FID on F1,F3 nuclei
%
% Note: decoupling with respect to any of the working nuclei is 
%       impossible in this pulse sequence.
%
% Andrew Porter
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=pansy_triple.m>

function fid=pansy_triple(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Timestep
timesteps=1./parameters.sweep;

% Initial state
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');

% Detection state
coil_h=state(spin_system,'L+',parameters.spins{1},'cheap');
coil_x1=state(spin_system,'L+',parameters.spins{2},'cheap');
coil_x2=state(spin_system,'L+',parameters.spins{3},'cheap');

% Get pulse operators
Hp=operator(spin_system,'L+',parameters.spins{1});
x1p=operator(spin_system,'L+',parameters.spins{2});
x2p=operator(spin_system,'L+',parameters.spins{3});
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; x1y=(x1p-x1p')/2i;
x2y=(x2p-x2p')/2i;

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
rho_stack_p=step(spin_system,Hy+x1y+x2y,rho_stack_p,pi/2);
rho_stack_m=step(spin_system,Hy+x1y+x2y,rho_stack_m,pi/2);

% Perform axial peak elimination
rho_stack=rho_stack_p-rho_stack_m;

% Run F2 evolution
fid.aa=evolution(spin_system,L,coil_h,rho_stack,timesteps(1),...
                 parameters.npoints(1)-1,'observable');
fid.ab=evolution(spin_system,L,coil_x1,rho_stack,timesteps(2),...
                 parameters.npoints(2)-1,'observable');
fid.ac=evolution(spin_system,L,coil_x2,rho_stack,timesteps(3),...
                 parameters.npoints(3)-1,'observable');             

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
elseif numel(parameters.sweep)~=3
    error('parameters array should have exactly three elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=3
    error('parameters.spins cell array should have exactly three elements.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.points variable.');
elseif numel(parameters.npoints)~=3
    error('parameters.npoints array should have exactly three elements.');
end
end

% Eternal war. We only dream of peace;
% May nothing stir us from the reverie.
% The hoary night, the icy cobalt bliss
% Of sleepy crows swaying in a tree.
%
% Joseph Brodsky, translated by IK

