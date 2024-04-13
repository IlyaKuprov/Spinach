% Cross-polarisation experiment in the rotating frame, followed by
% time-domain FID acquisition. The CP stage is preceded by wiping
% of the low-gamma spins and followed by FID acquisition with deco-
% upling of the high-gamma spins. Syntax:
%
%        fid=cp_acquire_soft(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.spins      - working spins, a cell array of
%                             strings with high-gamma spin fi-
%                             rst, and low-gamma spin last,
%                             for example {'1H','13C'}
%
%     parameters.hi_pwr     - nutation frequency of the exci-
%                             tation pulse on the high-gamma
%                             channel, Hz
%
%     parameters.cp_pwr     - nutation frequencies on the two
%                             channels during the CP contact
%                             time, a two-element vector, Hz
%
%     parameters.cp_dur     - duration of the contact time, s
%
%     parameters.rho0       - initial state, the state of the
%                             low-gamma spins will be wiped
%
%     parameters.coil       - detection state
%
%     H - Hamiltonian matrix, received from context function
%
%     R - relaxation superoperator, received from context function
%
%     K - kinetics superoperator, received from context function
%
% Output:
%
%     fid - signal detected on the coil state during 
%           system evolution
%
% marina.carravetta@soton.ac.uk
% p.t.williamson@soton.ac.uk
% guinevere.mathies@uni-konstanz.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cp_acquire_soft.m>

function fid=cp_acquire_soft(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Wipe the state of 13C
[~,rho]=decouple(spin_system,[],parameters.rho0,...
                                parameters.spins(2));

% Build and project 1H and 13C control operators
Hx=operator(spin_system,'Lx',parameters.spins{1});
Hy=operator(spin_system,'Ly',parameters.spins{1});
Cx=operator(spin_system,'Lx',parameters.spins{2});
Hx=kron(speye(parameters.spc_dim),Hx);
Hy=kron(speye(parameters.spc_dim),Hy);
Cx=kron(speye(parameters.spc_dim),Cx);

% Apply the 90-degree pulse on 1H along +X
rho=step(spin_system,L+2*pi*parameters.hi_pwr*Hx,...
                     rho,1/(4*parameters.hi_pwr));

% Run the CP contact time evolution: irradiation 
% of 1H along -Y, and of 13C along +X 
rho=evolution(spin_system,L-2*pi*parameters.cp_pwr(1)*Hy...
                           +2*pi*parameters.cp_pwr(2)*Cx,...
              [],rho,parameters.cp_dur,1,'final');

% Wipe the state of 1H and apply 1H decoupling
[L,rho]=decouple(spin_system,L,rho,parameters.spins(1));

% Run the acquisition
fid=evolution(spin_system,L,parameters.coil,rho,...
              1/parameters.sweep,parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'hi_pwr')||(parameters.hi_pwr<=0)
    error('high RF amplitude must be specified in parameters.hi_pwr variable.');
end
if ~isfield(parameters,'cp_pwr')
    error('CP RF amplitudes must be specified in parameters.cp_pwr variable.');
end
if ~isfield(parameters,'cp_dur')
    error('CP duration must be specified in parameters.cp_dur variable.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
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
end

% "Dealing with drunk law students has taught 
%  me to keep up with legislation."
%
% A bouncer at a night 
% club in Oxford

