% Constant-time COSY pulse sequence with analytical coherence selec-
% tion from Figure 3a of (https://doi.org/10.1002/chem.201406283).
% Syntax:
%
%             fid=ct_cosy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep        a two-element vector giving the 
%                            sweep widths in F1 and F2 
%
%    parameters.npoints      a two-element vector giving the 
%                            number of points in F1 and F2 
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
%    fid        - two-dimensional free induction decay
%
% mrw1g16@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ct_cosy.m>

function fid=ct_cosy(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Initial condition
if ~isfield(parameters,'rho0')
    parameters.rho0=state(spin_system,'Lz',parameters.spins{1},'cheap');
end

% Detection state
if ~isfield(parameters,'coil')
    parameters.coil=state(spin_system,'L+',parameters.spins{1},'cheap');
end

% Get the time grid for the CT period
CT=parameters.npoints(1)/parameters.sweep(1);
t1_grid=linspace(0,CT,parameters.npoints(1));

% Get the pulse operator
Hp=operator(spin_system,'L+',parameters.spins{1});

% Apply the first 90 degree pulse
rho=step(spin_system,(Hp+Hp')/2,parameters.rho0,pi/2);

% Select "+1" coherence
rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

% Preallocate rho stack
rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i);

% Loop over the value of t1
parfor n=1:parameters.npoints(1)
    
    % Run the first delay
    rho_current=evolution(spin_system,L,[],rho,CT/2-t1_grid(n)/2,1,'final');
    
    % Run the pi pulse
    rho_current=step(spin_system,(Hp+Hp')/2,rho_current,pi);
    
    % Run the second delay
    rho_current=evolution(spin_system,L,[],rho_current,CT/2+t1_grid(n)/2,1,'final');
    
    % Assign the stack element
    rho_stack(:,n)=rho_current;
    
end

% Final 90 degree pulse
rho_stack=step(spin_system,(Hp+Hp')/2,rho_stack,pi/2);

% Run the F2 evolution
fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.sweep(2),...
              parameters.npoints(2)-1,'observable');
          
% Flip t1 direction
fid=fliplr(fid);

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
    error('parameters.sweep array should have exactly two elements.');
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
end

% "Thanks to Professor Berry, who without none of this, 
%  would have been possible."
%
% Student feedback received by 
% Michael Berry on one of his
% modules

