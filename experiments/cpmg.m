% CPMG echo train with detection. Syntax:
%
%          fid=cpmg(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.nloops   - number of CPMG loops
%
%    parameters.timestep - time step
%
%    parameters.npoints  - number of steps per half-echo
%
% Outputs:
%
%    fid - free induction decay throughout the sequence
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cpmg.m>

function fid=cpmg(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Project the operator
parameters.pulse_op=kron(speye(parameters.spc_dim),parameters.pulse_op);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Excitation pulse
rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/2);

% Run a half-echo
traj=evolution(spin_system,L,[],rho,...
     parameters.timestep,parameters.npoints-1,'trajectory');

% Record the half-echo
fid=parameters.coil'*traj;

% CPMG loop
for n=1:parameters.nloops
    
    % Apply the pulse
    traj(:,end)=step(spin_system,parameters.pulse_op,traj(:,end),pi);
    
    % Run the echo
    traj=evolution(spin_system,L,[],traj(:,end),...
         parameters.timestep,2*parameters.npoints-1,'trajectory');
     
    % Record the echo
    fid=[fid parameters.coil'*traj]; %#ok<AGROW>
    
end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K) %#ok<INUSL>
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'pulse_op')
    error('pulse operator must be specified in parameters.pulse_op variable.');
end
end

% "Listen, man-cub," said the Bear, and his voice rumbled like thunder on a hot 
% night. "I have taught thee all the Law of the Jungle for all the peoples of 
% the jungle - except the Monkey-Folk who live in the trees. They have no law. 
% They are outcasts. They have no speech of their own, but use the stolen words 
% which they overhear when they listen, and peep, and wait up above in the bran-
% ches. Their way is not our way. They are without leaders. They have no remem-
% brance. They boast and chatter and pretend that they are a great people about
% to do great affairs in the Jungle, but the falling of a nut turns their minds
% to laughter and all is forgotten [...] and they desire, if they have any fixed
% desire, to be noticed by the Jungle People. But we do not notice them even
% when they throw nuts and filth on our heads."
%
% Rudyard Kiling, "The Jungle Book"

