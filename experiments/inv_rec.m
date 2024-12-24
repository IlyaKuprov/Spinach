% Inversion-recovery pulse sequence. Syntax:
%
%             fids=inv_rec(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep              spectrum sweep width, Hz
%
%    parameters.npoints            number of points in the FID
%
%    parameters.spins              nuclei on which the sequence runs,
%                                  specified as {'1H'}, {'13C'}, etc.
%
%    parameters.max_delay          longest relaxation delay
%
%    parameters.n_delays           number of relaxation delays to run
%
%    H     - Hamiltonian matrix, received from context function
%
%    R     - relaxation superoperator, received from context function
%
%    K     - kinetics superoperator, received from context function
%
% Outputs:
%
%    fids  - free induction decays for each delay starting from zero,
%            a matrix with individual FIDs in columns
%
% Note: the relaxation superoperator must be thermalised.
%
% Zak El-Machachi
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=inv_rec.m>

function fids=inv_rec(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get the thermodynamic equilibrium
H=hamiltonian(assume(spin_system,'labframe'),'left');
rho0=equilibrium(spin_system,H);

% Get the detection state
coil=state(spin_system,'L+',parameters.spins{1});

% Get the pulse operator
Lp=operator(spin_system,'L+',parameters.spins{1}); Ly=(Lp-Lp')/2i;

% Run a 180-degree pulse
rho=step(spin_system,Ly,rho0,pi);

% Run relaxation periods
rho_stack=evolution(spin_system,L,[],rho,...
                    parameters.max_delay/parameters.n_delays,...
                    parameters.n_delays,'trajectory');

% Run a 90-degree pulse
rho_stack=step(spin_system,Ly,rho_stack,pi/2);

% Run the detection period 
fids=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep,...
               parameters.npoints-1,'observable');
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
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if strcmp(spin_system.rlx.equilibrium,'zero')
    error('inter.equilibrium cannot be ''zero'' in this experiment.');
end
end

% Suddenly, everybody knows, and nothing's changed.
%
% Bruce Schneier, 
% about NSA surveillance

