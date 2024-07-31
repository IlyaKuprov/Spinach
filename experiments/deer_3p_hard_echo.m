% Samples the spin echo in the three-pulse DEER experiment to determine
% its precise location - in simulations, DEER spin echoes can be narrow
% and easy to miss. For high-spin electrons the echo may not be in the 
% expected place. Ideal hard pulses are used, each pulse only hits its
% specific electron or transition, depending on the pulse operator supp-
% lied. Syntax:
%
%          echo=deer_3p_hard_echo(spin_system,parameters,H,R,K)
%
% Parameters:
%
%      parameters.ex_prob   - probe pulse operator
%
%      parameters.ex_pump   - pump pulse operator
%
%      parameters.ta        - time between first and third pulse, s 
%
%      parameters.tb        - time between first and second pulse, s
%
%      parameters.tc        - time interval to sample around the 
%                             echo, s
%
%      parameters.rho0      - initial condition
%
%      parameters.nsteps    - number of sampling steps in the tc 
%                             interval
%
%      parameters.coil      - detection state
%
%      H  - Hamiltonian matrix, received from context function
%
%      R  - relaxation superoperator, received from context function
%
%      K  - kinetics superoperator, received from context function
%
% Outputs:
%
%      echo - the signal detected during the parameters.tc interval
%             with parameters.nsteps points in it
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
% nurit.manukovsky@weizmann.ac.il
% daniella.goldfarb@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=deer_3p_hard_echo.m>

function echo=deer_3p_hard_echo(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K)

% Compose Liouvillian
L=H+1i*R+1i*K;

% First pulse
rho=step(spin_system,parameters.ex_prob,parameters.rho0,pi/2);

% Evolution
rho=evolution(spin_system,L,[],rho,parameters.tb,1,'final');

% Second pulse
rho=step(spin_system,parameters.ex_pump,rho,pi);
     
% Evolution
rho=evolution(spin_system,L,[],rho,parameters.ta-parameters.tb,1,'final');

% Third pulse
rho=step(spin_system,parameters.ex_prob,rho,pi);

% Evolution
rho=evolution(spin_system,L,[],rho,parameters.ta-0.5*parameters.tc,1,'final');

% Evolution
echo=evolution(spin_system,L,parameters.coil,rho,...
               parameters.tc/parameters.nsteps,parameters.nsteps,'observable');

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
    error('number of steps must be specificed in parameters.steps variable.');
end
if (~isnumeric(parameters.nsteps))||(~isreal(parameters.nsteps))||...
   (numel(parameters.nsteps)~=1)||(parameters.nsteps<1)||(mod(parameters.nsteps,1)~=0)
    error('parameters.nsteps must be a positive integer.');
end
if ~isfield(parameters,'ta')
    error('time between first and third pulse must be specificed in parameters.ta variable.');
end
if (~isnumeric(parameters.ta))||(~isreal(parameters.ta))||...
   (numel(parameters.ta)~=1)||(parameters.ta<=0)
    error('parameters.ta must be a positive number.');
end
if ~isfield(parameters,'tb')
    error('time between first and second pulse must be specificed in parameters.tb variable.');
end
if (~isnumeric(parameters.tb))||(~isreal(parameters.tb))||...
   (numel(parameters.tb)~=1)||(parameters.tb<=0)
    error('parameters.tb must be a positive number.');
end
if ~isfield(parameters,'tc')
    error('time to sample around the echo must be specificed in parameters.tc variable.');
end
if (~isnumeric(parameters.tc))||(~isreal(parameters.tc))||...
   (numel(parameters.tc)~=1)||(parameters.tc<=0)
    error('parameters.tc must be a positive number.');
end
if ~isfield(parameters,'rho0')
    error('the initial state must be provided in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('the detection state must be provided in parameters.coil variable.');
end
if ~isfield(parameters,'ex_prob')
    error('probe pulse operator must be provided in parameters.ex_prob variable.');
end
if ~isfield(parameters,'ex_pump')
    error('pump pulse operator must be provided in parameters.ex_pump variable.');
end
if (0.5*parameters.tc>parameters.ta)||(parameters.tb>parameters.ta)
    error('incorrect timing settings.');
end
end

% History is a set of lies agreed upon.
%
% Napoleon Bonaparte

