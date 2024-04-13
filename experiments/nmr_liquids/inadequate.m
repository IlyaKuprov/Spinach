% INADEQUATE pulse sequence. Removes uncoupled carbons from the spectra
% via a double quantum filter, allowing only double quantum coherences
% states too be detected. At natural abundance 13C, this produces only
% 13C pair subspectra. Syntax:
%
%              fid=inadequate(spin_system,parameters,H,R,K)
%
% Parameters: 
% 
%    parameters.sweep              sweep width in Hz
%
%    parameters.npoints            number of points in the fid
%
%    parameters.spins              active nuclei, e.g. {'13C'}
%
%    parameters.J                  working J-coupling in Hz
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Output:
%
%    fid  - free induction decay
%
% Note: use dilute.m to generate carbon pair isotopomers.
% 
% Bud Macaulay, Ilya Kuprov
%
% <https://spindynamics.org/wiki/index.php?title=inadequate.m>

function fid=inadequate(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Decoupling
L=decouple(spin_system,L,[],parameters.decouple);

% Sequence timing
timestep=1./parameters.sweep;
tau=abs(1/(4*parameters.J));

% Initial and detection states
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Pulse operators
Cp=operator(spin_system,'L+',parameters.spins{1});
Cx=(Cp+Cp')/2; Cy=(Cp-Cp')/2i;

% Pulse 90x
rho=step(spin_system,Cx,rho,pi/2);

% J-coupling evolution
rho=step(spin_system,L,rho,tau);

% Pulse 180y
rho=step(spin_system,Cy,rho,pi);

% J-coupling evolution
rho=step(spin_system,L,rho,tau);

% Pulse 90x 
rho=step(spin_system,Cx,rho,pi/2);

% Select DQC
rho=coherence(spin_system,rho,{{parameters.spins{1},[2 -2]}});

% Pulse on 90x
rho=step(spin_system,Cx,rho,pi/2);

% Detection
fid=evolution(spin_system,L,coil,rho,timestep(1),...
              parameters.npoints(1)-1,'observable');

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
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=1
    error('parameters.npoints array should have exactly one element.');
end
if ~isfield(parameters,'J')
    error('scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('parameters.J array should have exactly one element.');
end
end

% When children ask "why do beavers have big teeth", adults usually say
% "to help with felling trees and building dams." However, the correct 
% answer is "because beavers with big teeth survived, and beavers with
% small teeth did not". That's how evolution works - we see the species
% that had survived, and their features are those that helped them sur-
% vive. Nature does not assist animals, it only spares those whose mu-
% tations were more fortuitous. Thus, any question about animals may be
% answered quickly: "others died". "Why do rabbits become white for win-
% ter?" - "Wolves ate the grey ones." - "Why do zebras have stripes?" -
% "Lions killed the ones that had no stripes." - "Why do elephants have
% big ears?" - "Elephants with small ears died from heat stroke" - "Why
% do capybaras..." - because the rest perished from disease, got eaten,
% died of hunger, fell to death, drowned... And if you don't stop pick-
% ing your nose, your brother will be the one who survives.
%
% A Soviet biology textbook
% for primary schools

