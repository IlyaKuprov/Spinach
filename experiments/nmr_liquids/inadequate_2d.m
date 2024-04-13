% 2D INADEQUATE pulse sequence. Syntax:
%
%         fid=inadequate_2d(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.sweep              [F1 F2] sweep widths, Hz
%
%     parameters.npoints            [F1 F2] numbers of points
%
%     parameters.spins              {F1 F2} nuclei (e.g. '13C','1H')
%
%     parameters.J                  working scalar coupling, Hz
%
%     H  - Hamiltonian matrix, received from context function
%
%     R  - relaxation superoperator, received from context function
%
%     K  - kinetics superoperator, received from context function
%
% Outputs:
%
%     fid.cos,fid.sin -  two components of the States signal
%
% Notes: use dilute.m to generate carbon pair isotopomers.
%
% Theresa Hune
% Christian Griesinger
%
% <https://spindynamics.org/wiki/index.php?title=inadequate_2d.m>

function fid=inadequate_2d(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get timing parameters
tau=abs(1/(4*parameters.J));
timestep=1./parameters.sweep;

% Initial state along Z
rho=state(spin_system,'Lz',parameters.spins{1});

% Quadrature detection state
coil=state(spin_system,'L+',parameters.spins{1});

% Pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;
Lxy=cos(pi/4)*Lx+sin(pi/4)*Ly;

% 90x pulse
rho_cos=step(spin_system,Lx,rho,pi/2);
rho_sin=step(spin_system,Lxy,rho,pi/2);

% J-coupling evolution during tau
rho_cos=evolution(spin_system,L,[],rho_cos,tau,1,'final');
rho_sin=evolution(spin_system,L,[],rho_sin,tau,1,'final');

% inversion pulse
rho_cos=step(spin_system,Lx,rho_cos,pi);
rho_sin=step(spin_system,Lxy,rho_sin,pi);

% Second tau evolution
rho_cos=evolution(spin_system,L,[],rho_cos,tau,1,'final');
rho_sin=evolution(spin_system,L,[],rho_sin,tau,1,'final');

% 90x pulse, creation of MQC
rho_cos=step(spin_system,Lx,rho_cos,pi/2);
rho_sin=step(spin_system,Lxy,rho_sin,pi/2);

% Apply the double-quantum filter
rho_cos=coherence(spin_system,rho_cos,{{parameters.spins{1},[+2,-2]}});
rho_sin=coherence(spin_system,rho_sin,{{parameters.spins{1},[+2,-2]}});

% t1 evolution
rho_stack_cos=evolution(spin_system,L,[],rho_cos,timestep(1),parameters.npoints(1)-1,'trajectory');
rho_stack_sin=evolution(spin_system,L,[],rho_sin,timestep(1),parameters.npoints(1)-1,'trajectory');

% 90x pulse, creation of observable magnetization
rho_stack_cos=step(spin_system,Lx,rho_stack_cos,pi/2);
rho_stack_sin=step(spin_system,Lx,rho_stack_sin,pi/2);

% Run the F2 evolution
fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable');
fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable');

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
if ~isfield(parameters,'J')
    error('scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('parameters.J array should have exactly one element.');
end
end

% Experiments are the only means of knowledge at 
% our disposal. The rest is poetry, imagination.
%
% Max Planck

