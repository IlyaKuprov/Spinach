% ESEEM pulse sequence with ideal hard pulses. Syntax:
%
%             fid=eseem(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.npoints            number of points to be computed
%
%    parameters.timestep           simulation time step, seconds
%
%    parameters.rho0               initial state
%
%    parameters.coil               detection state
%
%    parameters.screen             optional screen state (must be 
%                                  the Hermitian conjugate of the
%                                  detection state)
%
%    parameters.pulse_op           pulse operator
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid  - time domaing signal whose Fourier transform is the 
%           ESEEM spectrum
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=eseem.m>

function fid=eseem(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K)

% Compose Liouvillian
L=H+1i*R+1i*K;

% Set the defaults
if ~isfield(parameters,'screen'), parameters.screen=[]; end

% First pulse
rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/2);

% Spin echo
rho_stack=evolution(spin_system,L,[],rho,parameters.timestep/2,...
                   (parameters.npoints-1),'trajectory',parameters.screen);
% Second pulse
rho_stack=step(spin_system,parameters.pulse_op,rho_stack,pi);

% Spin echo
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep/2,...
                   (parameters.npoints-1),'refocus',parameters.coil);
% Detect
fid=transpose(full(parameters.coil'*rho_stack));

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if numel(parameters.npoints)~=1
    error('parameters.npoints array should have exactly one element.');
end
if ~isfield(parameters,'timestep')
    error('time step should be specified in parameters.timestep variable.');
end
if numel(parameters.timestep)~=1
    error('parameters.timestep array should have exactly one element.');
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
if isfield(parameters,'screen')&&(size(parameters.screen,1)~=size(parameters.rho0,1))
    error('if parameters.screen is specified, then it must have the same number of rows as parameters.rho0.');
end
if ~isfield(parameters,'pulse_op')
    error('pulse operator must be specified in parameters.pulse_op variable.');
end
if any(size(parameters.pulse_op)~=size(H))
    error('parameters.pulse_op must have the same dimension as the Hamiltonian.');
end
end

% In my life, I have got used to being called four-letter names.
% Going down to three is a completely new experience which I will
% hopefully enjoy.
%
% Sir Andre Geim, after receiving his knighthood in 2011.

