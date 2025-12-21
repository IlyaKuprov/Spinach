% DEPTQ pulse sequence from Figure 1 of the paper by Burger and Bigler:
% (https://doi.org/10.1006/jmre.1998.1595). Syntax: 
% 
%                  fid=deptq(spin_system,parameters,H,R,K)
%
% Parameters: 
% 
%    parameters.sweep              [F1] Sweep width in Hz
%
%    parameters.npoints            [F1] number of points
%
%    parameters.spins              {F1,F2} nuclei, e.g. {'13C','1H'}
%
%    parameters.J                  working J-coupling in Hz
% 
%    parameters.beta               the angle used in the selection
%                                  pulse, radians
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid  - free induction decay
%
% Note: use dilute.m to generate carbon isotopomers.
%
% Note: the sequence differs from dept.m in that quaternary carbons
%       do appear.
%
% Andrew Porter, Ilya Kuprov
%
% <https://spindynamics.org/wiki/index.php?title=deptq.m>

function fid=deptq(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K)

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get timing parameters
tau_j=abs(1/(2*parameters.J)); 
timestep=1/parameters.sweep;

% Isotropic thermal equilibrium
rho=equilibrium(spin_system);

% Detection state - Cy
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Pulse operators
Cp=operator(spin_system,'L+',parameters.spins{1}); 
Cx=(Cp+Cp')/2; Cy=(Cp-Cp')/2i;
Hp=operator(spin_system,'L+',parameters.spins{2}); 
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;

% Carbon pulse
rho=step(spin_system,Cy,rho,-pi/2);

% J-coupling evolution
rho=step(spin_system,L,rho,tau_j);

% Proton and carbon pulses
rho=step(spin_system,Hx,rho,+pi/2);
rho_a=step(spin_system,Cx,rho,+pi);
rho_b=step(spin_system,Cx,rho,-pi);
rho_c=step(spin_system,Cy,rho,+pi);
rho_d=step(spin_system,Cy,rho,-pi);

% J-coupling evolution
rho_a=step(spin_system,L,rho_a,tau_j);
rho_b=step(spin_system,L,rho_b,tau_j);
rho_c=step(spin_system,L,rho_c,tau_j);
rho_d=step(spin_system,L,rho_d,tau_j);

% Proton and carbon pulses
rho_a=+step(spin_system,Hx,rho_a,+pi)...
      -step(spin_system,Hy,rho_c,+pi); 
rho_b=+step(spin_system,Hx,rho_b,-pi)...
      -step(spin_system,Hy,rho_d,-pi); 
rho_a=+step(spin_system,Cx,rho_a,+pi/2);
rho_b=+step(spin_system,Cx,rho_b,+pi/2);

% J-coupling evolution
rho_a=step(spin_system,L,rho_a,tau_j);
rho_b=step(spin_system,L,rho_b,tau_j);

% Proton and carbon pulses
rho=step(spin_system,Cx,rho_a,+pi)+...
    step(spin_system,Cx,rho_b,-pi);
rho=step(spin_system,Hy,rho,parameters.beta);

% J-coupling evolution
rho=step(spin_system,L,rho,tau_j);

% Detection pulse
rho=step(spin_system,Cx,rho,-pi/2);

% Proton decoupling
[L,rho]=decouple(spin_system,L,rho,parameters.spins(2));

% Phase cycled detection
fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('dept: sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('dept: parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('dept: working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('dept: parameters.spins cell array should have exactly two elements.');
end
if ~isfield(parameters,'npoints')
    error('deptt: number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=1
    error('dept: parameters.npoints array should have exactly one element.');
end
if ~isfield(parameters,'J')
    error('dept: scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('dept: parameters.J array should have exactly one element.');
end
if ~isfield(parameters,'beta')
    error('dept: Beta angle should be specified in parameters.beta variable.');
elseif numel(parameters.beta)~=1
    error('dept: parameters.beta array should have exactly one element.');
end
end

% After publication of "The Hobbit" in 1937, the grammar pedants 
% of the time slammed Tolkien for using "dwarves" as the plural
% for "dwarf". They were pointing out that, in the Oxford English
% Dictionary, the only correct form is "dwarfs". Tolkien ignored
% the pedants. He was the Editor of the Oxford English Dictionary.

