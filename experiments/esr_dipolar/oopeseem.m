% Out-of-phase ESEEM pulse sequence with the first pulse set to pi/4 to
% probe two-electron correlations in the initial condition. Syntax:
%
%             fid=eseem(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.npoints            number of time points to be 
%                                  computed
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
%    parameters.pulse_op           pulse operator A, the propagators
%                                  for the pulses will be exp(-i*A*pi)
%                                  and exp(-i*A*pi/4)
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid   -  OOP-ESEEM time trace
%
% Note: the sequence uses ideal pulses, replace with shaped_pulse_af()
%       to have soft pulses instead.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=oopeseem.m>

function fid=oopeseem(spin_system,parameters,H,R,K)

% Move into adjoint representation if needed
if strcmp(spin_system.bas.formalism,'zeeman-hilb')

    % Inform the user
    report(spin_system,'projecting zeeman-hilb simulation into Liouville space...');

    % Project into Liouville space
    H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm');

    % Stretch the state parameters
    if isfield(parameters,'rho0')
        parameters.rho0=reshape(parameters.rho0,size(spin_system.bas.basis,1)^2,[]);
    end
    if isfield(parameters,'coil')
        parameters.coil=reshape(parameters.coil,size(spin_system.bas.basis,1)^2,[]);
    end

    % Stretch the screen state
    if isfield(parameters,'screen')
        parameters.screen=reshape(parameters.screen,size(spin_system.bas.basis,1)^2,[]);
    end

    % Project the pulse operator
    if isfield(parameters,'pulse_op')
        parameters.pulse_op=hilb2liouv(parameters.pulse_op,'comm');
    end

    % Rebuild the basis index table for the Liouville space
    zbas=spin_system.bas.basis; hdim=size(zbas,1);
    spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))];

    % Discard Hilbert space symmetry data
    if isfield(spin_system.bas,'irrep')
        spin_system.bas=rmfield(spin_system.bas,'irrep');
        report(spin_system,'Hilbert space irreps discarded, proceeding without symmetry.');
    end

    % Update the formalism setting
    spin_system.bas.formalism='zeeman-liouv';

end

% Check consistency
grumble(spin_system,parameters,H,R,K)

% Compose Liouvillian
L=H+1i*R+1i*K;

% Set the defaults
if ~isfield(parameters,'screen'), parameters.screen=[]; end

% First pulse
rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/4);

% Spin echo
rho_stack=evolution(spin_system,L,[],rho,parameters.timestep/2,...
                   (parameters.npoints-1),'trajectory',parameters.screen);
% Second pulse
rho_stack=step(spin_system,parameters.pulse_op,rho_stack,pi);

% Spin echo
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep/2,...
                   (parameters.npoints-1),'refocus',parameters.coil);
% Detect
fid=transpose(parameters.coil'*rho_stack);

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
    error('parameters.screen must have the same number of rows as parameters.rho0.');
end
if ~isfield(parameters,'pulse_op')
    error('pulse operator must be specified in parameters.pulse_op variable.');
end
if any(size(parameters.pulse_op)~=size(H))
    error('parameters.pulse_op must have the same dimension as the Hamiltonian.');
end
end

% If you want to cut your own throat, don't come
% to me for a bandage.
%
% Margaret Thatcher

