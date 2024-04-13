% INEPT pulse sequence. Syntax:
% 
%           fid=inept(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep           [F1] Sweep width in Hz
%
%    parameters.npoints         [F1] number of points
%
%    parameters.spins           {F1 F2} working nuclei,
%                               e.g. {'15N','1H'}
%
%    parameters.J               working scalar coupling 
%                               in Hz
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid   - free induction decay
%
% Note: use dilute() to generate carbon isotopomers.
%
% Andrew Porter, Ilya Kuprov
%
% <https://spindynamics.org/wiki/index.php?title=inept.m>

function fid=inept(spin_system,parameters,H,R,K)

% consistency check
grumble(spin_system,parameters,H,R,K)

% compose Liouvillian
L=H+1i*R+1i*K;

% Get timing parameters
tau=abs(1/(4*parameters.J));
timestep=1/parameters.sweep;

% Get accurate thermal equilibrium
H0=hamiltonian(assume(spin_system,'labframe'),'left');
rho=equilibrium(spin_system,H0);

% Detection state
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Pulse operators
Cp=operator(spin_system,'L+',parameters.spins{1}); 
Cx=(Cp+Cp')/2; Cy=(Cp-Cp')/2i;
Hp=operator(spin_system,'L+',parameters.spins{2}); 
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;

% 90x pulse on H
rho=step(spin_system,Hx,rho,pi/2);

% tau evolution
rho=evolution(spin_system,L,[],rho,tau,1,'final');

% Two inversion pulses
rho=step(spin_system,Cy+Hy,rho,pi);

% Second tau evolution
rho=evolution(spin_system,L,[],rho,tau,1,'final');

% 90x pulse on C
rho=step(spin_system,Cx,rho,pi/2);

% Split phase 90y pulses on H
rho_pos=step(spin_system,Hy,rho,+pi/2);
rho_neg=step(spin_system,Hy,rho,-pi/2);

% Phase cycle
rho=(rho_pos-rho_neg)/2;

% Detection
fid=evolution(spin_system,L,coil,rho,timestep,...
              parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('inept: sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('inept: parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('inept: working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('inept: parameters.spins cell array should have exactly two elements.');
end
if ~isfield(parameters,'npoints')
    error('inept: number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=1
    error('inept: parameters.npoints array should have exactly one element.');
end
if ~isfield(parameters,'J')
    error('inept: scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('inept: parameters.J array should have exactly one element.');
end
end

% Two dogs are walking along a street. They are passed by a third dog
% driving a lorry load of logs. One turns to the other and says "This
% guy started fetching a stick and built up the business from there."
%
% An American joke

