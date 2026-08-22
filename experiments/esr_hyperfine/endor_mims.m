% Mims ENDOR pulse sequence with ideal hard pulses. Syntax:
%
%          fid=endor_mims(spin_system,parameters,H,R,K)
%
% Parameters:
%
%   parameters.sweep           nuclear frequency sweep width, Hz
%
%   parameters.npoints         number of fid points to be computed
%
%   parameters.tau             stimulated echo time, seconds
%
%   H  - Hamiltonian matrix, received from context function
%
%   R  - relaxation superoperator, received from context function
%
%   K  - kinetics superoperator, received from context function
%
% Outputs:
%
%   fid  - free induction decay whose Fourier transform is the 
%          Mims ENDOR signal
%
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=endor_mims.m>

function fid=endor_mims(spin_system,parameters,H,R,K)

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

    % Rebuild the basis index table for the Liouville space
    zbas=spin_system.bas.basis; hdim=size(zbas,1);
    spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))];

    % Update the formalism setting
    spin_system.bas.formalism='zeeman-liouv';

end

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Initial state
rho=state(spin_system,'Lz','electrons');

% Detection state
coil=state(spin_system,'L+','electrons');

% Pulse operators
Ep=operator(spin_system,'L+','electrons'); Ey=(Ep-Ep')/2i;
Np=operator(spin_system,'L+','nuclei'); Ny=(Np-Np')/2i;
        
% Apply the initial pulses
rho=step(spin_system,Ey,rho,pi/2);
rho=evolution(spin_system,L,[],rho,parameters.tau,1,'final');
rho=step(spin_system,Ey,rho,pi/2);

% Analytical coherence selection
rho=coherence(spin_system,rho,{{'electrons',0}});

% Apply pulses on nuclear spins
rho=step(spin_system,Ny,rho,pi/2);

% Analytical coherence selection
rho=coherence(spin_system,rho,{{'nuclei',[-1 1]}});

% Indirect dimension evolution
rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep,...
                    parameters.npoints-1,'trajectory');
rho_stack=+step(spin_system,+Ny,rho_stack,pi/2) ...
          -step(spin_system,-Ny,rho_stack,pi/2);

% Apply pulse on electron spin to refocus
rho_stack=step(spin_system,Ey,rho_stack,pi/2);
rho_stack=evolution(spin_system,L,[],rho_stack,...
                    parameters.tau,1,'final',coil);

% Detect
fid=full(transpose(coil'*rho_stack));

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
if ~isfield(parameters,'sweep')
    error('endor_mims: sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('endor_mims: parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('endor_mims: number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=1
    error('endor_mims: parameters.npoints array should have exactly one element.');
end
if ~isfield(parameters,'tau')
    error('endor_mims: echo time should be specified in parameters.tau variable.');
elseif (~isnumeric(parameters.tau))||(~isreal(parameters.tau))||...
       (numel(parameters.tau)~=1)||(~isfinite(parameters.tau))||...
       (parameters.tau<=0)
    error('endor_mims: parameters.tau must be a finite positive real scalar.');
end
end

% Reason is not automatic. Those who deny it cannot be conquered
% by it. Do not count on them. Leave them alone.
%
% Ayn Rand

