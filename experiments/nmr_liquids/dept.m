% DEPT pulse sequence from the schematic in the paper by Doddrell, Pegg,
% and Bendall (https://doi.org/10.1016/0022-2364(82)90286-4). Syntax: 
% 
%                  fid=dept(spin_system,parameters,H,R,K)
%
% Parameters: 
% 
%    parameters.sweep          [F1] Sweep width in Hz
%
%    parameters.npoints        [F1] number of points
%
%    parameters.spins          {F1,F2} nuclei, e.g. {'13C','1H'}
%
%    parameters.J              working J-coupling in Hz
% 
%    parameters.beta           the angle used in the selection
%                              pulse, radians
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
% Note: DEPT135 yields spectra with CH and CH3 signals in opposite phase
%       to CH2 signals; DEPT90 yields spectra with only CH signals; DEPT45
%       yields spectra with positive CH, CH2, and CH3 signals; quaternary
%       carbons do not appear.
%
% Note: use dilute.m to generate carbon isotopomers.
%
% a.porter@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dept.m>

function fid=dept(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get timing parameters
tau_j=abs(1/(2*parameters.J)); 
timestep=1/parameters.sweep;

% Get accurate thermal equilibrium
H0=hamiltonian(assume(spin_system,'labframe'),'left');
rho=equilibrium(spin_system,H0);

% Detection state
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Pulse operators
Cp=operator(spin_system,'L+',parameters.spins{1});
Hp=operator(spin_system,'L+',parameters.spins{2}); 
Cx=(Cp+Cp')/2; Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;

% Proton pulse
rho=step(spin_system,Hx,rho,pi/2);

% J-coupling evolution
rho=step(spin_system,L,rho,tau_j);

% Proton and carbon pulses
rho=step(spin_system,Cx,rho,pi/2);
rho_a=step(spin_system,Hx,rho,+pi)-step(spin_system,Hy,rho,+pi);
rho_b=step(spin_system,Hx,rho,-pi)-step(spin_system,Hy,rho,-pi);

% Second tau evolution
rho_a=step(spin_system,L,rho_a,tau_j);
rho_b=step(spin_system,L,rho_b,tau_j);

% Proton and carbon pulses
rho=step(spin_system,Cx,rho_a,+pi)+...
    step(spin_system,Cx,rho_b,-pi);
rho=step(spin_system,Hy,rho,parameters.beta);

% Third tau evolution
rho=step(spin_system,L,rho,tau_j);

% Decoupling
[L,rho]=decouple(spin_system,L,rho,parameters.spins(2));

% Detection
fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable');

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
    error('dept: beta angle should be specified in parameters.beta variable.');
elseif numel(parameters.beta)~=1
    error('dept: parameters.beta array should have exactly one element.');
end
end

% IK was recently asked why scientists hate working with other people's
% code. After a lot of thinking on how to cast the nightmare into terms 
% that a lay person would understand, here's a good example.
%
% Say, you are asked to complete another builder's unfinished project -
% a research lab on an uninhabited island. You get to the site and see,
% apart from the unfinished building, a big (building-sized) fan and a
% hot air balloon next to it. In the building, there is a room stuffed
% full of broomsticks. You scratch your head, remove all of that weird
% nonsense and finish the building. 
%
% Five minutes after the building is handed over to scientists, they run
% out and scream "POISON GAS LEAK!!!". "What the hell", you think, "it
% must all just work!" You phone the previous builder: "Jim, we have a
% poson gas leak. What's the problem?" - "No idea, everything should be
% fine. Did you change anything in the project?" - "Not much, we threw
% the broomsticks away..." - "The broomsticks were there to hold up the 
% ceiling!" - "Dafuq?!" - "They were holding the ceiling, I am telling
% you. The gas cistern is directly above them. Very heavy, so we had to
% reinforce with broomsticks." - "You should have left a bloody note at
% least! OK, we have a gas leak, what do we do?" - "Turn on the fan, it
% will blow the gas away from the island." - "I had dismounted that fan
% the moment I got here!" - "Why?" - "And why did you build a 120-tonne
% fan? Could you not just put a box of bloody gas masks?" - "A box of
% gas masks you have to look for, and I had that fan left over from a
% previous project." - "Jim, I had removed your fan, we are suffocating
% here!" - "What the hell are you still doing there? Climb into the hot
% air balloon and get the hell out!"

