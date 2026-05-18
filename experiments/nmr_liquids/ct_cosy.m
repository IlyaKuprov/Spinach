% Constant-time COSY pulse sequence with analytical coherence selec-
% tion.
%
% References:
%
%             https://doi.org/10.1006/jmra.1994.1095
%             https://doi.org/10.1080/00387010009350054
%             https://doi.org/10.1002/chem.201406283
%
% The F1 trace is flipped at the end to preserve the conventional
% indirect-dimension sign for the reversed-delay implementation.
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
%    parameters.angle        final pulse angle in radians, defaults
%                            to pi/2
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

% Set default final pulse angle
if ~isfield(parameters,'angle')
    parameters.angle=pi/2;
end

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
t1_grid=(0:(parameters.npoints(1)-1))/parameters.sweep(1);
CT=t1_grid(end);

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

% Final pulse
rho_stack=step(spin_system,(Hp+Hp')/2,rho_stack,parameters.angle);

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
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       any(~isfinite(parameters.sweep))||any(parameters.sweep<=0)
    error('parameters.sweep must contain two positive real numbers.');
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
    error('final pulse angle should be specified in parameters.angle variable.');
elseif numel(parameters.angle)~=1
    error('parameters.angle array should have exactly one element.');
elseif (~isnumeric(parameters.angle))||(~isreal(parameters.angle))||...
       (~isfinite(parameters.angle))
    error('parameters.angle must be a finite real scalar.');
end
end

% "Thanks to Professor Berry, who without none of this, 
%  would have been possible."
%
% Student feedback received by 
% Michael Berry on one of his
% modules

