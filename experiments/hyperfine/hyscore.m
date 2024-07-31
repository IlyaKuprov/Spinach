% HYSCORE experiment, implemented as described in Szosenfogel and
% Goldfarb (http://dx.doi.org/10.1080/00268979809483260). Syntax:
%
%             fid=hyscore(spin_system,parameters,H,R,K)
%
% The following parameters are currently accepted:
%
%    parameters.nsteps          number of points to be computed 
%                               in each dimension
%
%    parameters.sweep           sweep width, Hz
%
%    parameters.tau             tau delay, seconds
%
%    parameters.rho0            initial state
%
%    parameters.coil            detection state
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid  - two-dimensional free induction decay that Fourier
%           transforms into a HYSCORE spectrum
%
% Note: the sequence uses ideal pulses, replace with shaped_pulse_af()
%       to have soft pulses instead.
% 
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=hyscore.m>

function fid=hyscore(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get the pulse operators
Lp=operator(spin_system,'L+','E');
Lm=operator(spin_system,'L-','E');
Lx=(Lp+Lm)/2;

% Apply the first pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);

% Run the tau evolution
rho=evolution(spin_system,L,[],rho,parameters.tau,1,'final');

% Apply the second pulse
rho=step(spin_system,Lx,rho,pi/2);

% Apply coherence filter
rho=coherence(spin_system,rho,{{'E',0}});

% Run the indirect dimension evolution
rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep,...
                    parameters.nsteps(1)-1,'trajectory');

% Apply the third pulse
rho_stack=step(spin_system,Lx,rho_stack,pi);

% Propagate coil state backwards in time
coil=evolution(spin_system,L,[],parameters.coil,-parameters.tau,1,'final');

% Apply a backwards pulse on the coil
coil=step(spin_system,-Lx,coil,pi/2);

% Detect on new coil state in the direct dimension
fid=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep,...
              parameters.nsteps(2)-1,'observable');

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'nsteps')
    error('number of time steps should be specified in parameters.nsteps variable.');
end
if numel(parameters.nsteps)~=2
    error('parameters.nsteps array should have two elements.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if numel(parameters.sweep)~=1
    error('parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if size(parameters.rho0,1)~=size(H,2)
    error('dimensions of parameters.rho0 and Hamiltonian must be consistent.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if size(parameters.coil,1)~=size(H,2)
    error('dimensions of parameters.coil and Hamiltonian must be consistent.');
end
if ~isfield(parameters,'tau')
    error('tau delay should be specified in parameters.tau variable.');
end
if numel(parameters.tau)~=1
    error('parameters.tau array should have exactly one element.');
end
end

% The best summary IK has ever heard of the checks and regulations at UK
% universities is given by an old Russian joke: "A hitchhiker was picked
% up by a car on one of the roads in Siberia. After a few hundred miles,
% the driver was stopped by the police and slapped with an insignificant
% spot fine. The driver returned to the car, visibly irritated, but also
% entertained. Asked by the hitchhiker about what just happened, the dri-
% ver explained: 'What a bloody idiot. The driving license has long expi-
% red, the photo in it is not mine, the car is stolen, there is a packet
% of MDMA, an AK-47 and a dead body in the trunk -- and you know what he
% just fined me for? Not wearing the seat belt.'"

