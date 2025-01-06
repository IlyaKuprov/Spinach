% Magnetic field scan steady-state DNP experiment. Returns the 
% steady-state population of the user-specified state as a fun-
% ction of magnetic field. Syntax:
%
%      dnp=dnp_field_scan(spin_system,parameters,H,R,K)
%
% Parameters:
%
%   parameters.mw_pwr      -   microwave power, Hz
%
%   parameters.mw_frq      -   microwave frequency offset from 
%                              the free electron frequency at 
%                              the reference B0 field, Hz
%
%   parameters.fields      -   a vector of magnetic field off-
%                              sets from the reference B0 field,
%                              Tesla
%
%   parameters.rho0        -   equilibrium state at the reference
%                              B0 field
%
%   parameters.coil        -   coil state vector or a horizon-
%                              tal stack thereof
%
%   parameters.mw_oper     -   microwave irradiation operator
%
%   parameters.ez_oper     -   Lz operator on the electrons
%
%   parameters.method      -   'backslash' to use Matlab's
%                              linear equation solver, 'gmres'
%                              to use ILU preconditioned GMRES
%
%   H - Hamiltonian matrix, received from context function
%
%   R - relaxation superoperator, received from context function
%
%   K - kinetics superoperator, received from context function
%
% Output:
%
%    dnp    -  an array of steady state expectation values for
%              the states specified in parameters.coil at each
%              of the fields supplied
% 
% Note: the relaxation superoperator should NOT be thermalized
%       for this type of calculation.
%
% Note: thermal equilibrium state and relaxation superoperator are
%       assumed to be the same at all fields in the sweep - DO NOT
%       USE with broad magnetic field sweep experiments.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=dnp_field_scan.m>

function dnp=dnp_field_scan(spin_system,parameters,H,R,K) 

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Damp unit state and check
R(1,1)=-mean(abs(diag(R))); Rc=condest(R);
if Rc>1e9, error('R must be non-singular.'); end
if ismember('redfield',spin_system.rlx.theories)&&...
           (Rc>1/spin_system.tols.rlx_integration)
    error('R acuracy too low, reduce spin_system.tols.rlx_integration');
end

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Preallocate the answer
dnp=zeros([numel(parameters.fields) size(parameters.coil,2)],'like',1i);

% Add microwave terms to the Liouvillian
L=L+2*pi*parameters.mw_pwr*parameters.mw_oper;

% Add offset terms to the Liouvillian
L=L-2*pi*parameters.mw_frq*parameters.ez_oper;

% Precompute the right hand side
b=R*parameters.rho0;

% Loop over the fields
parfor n=1:numel(parameters.fields)
    
    % Calculate electron offset frequency
    omega=-parameters.fields(n)*spin('E'); %#ok<PFBNS>
    
    % Add electron offset terms to the Liouvillian
    L_current=L+omega*parameters.ez_oper;
    
    % Get the steady state DNP
    if strcmp(parameters.method,'gmres')
        
        % Make sure L is sparse;
        L_current=sparse(L_current);
        
        % Update the preconditioner
        [M1,M2]=ilu(-1i*L_current,struct('type','crout','droptol',1e-6));
        
        % Run using preconditioned GMRES
        dnp(n,:)=parameters.coil'*gmres(-1i*L_current,b,10,1e-10,numel(parameters.rho0),M1,M2);
        
    elseif strcmp(parameters.method,'backslash')
        
        % Run using Matlab's backslash
        dnp(n,:)=parameters.coil'*(-1i*L_current\b);
        
    else
        
        % Complain and bomb out
        error('unknown solver.');
        
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'mw_pwr')
    error('microwave power should be specified in parameters.mw_pwr variable.');
elseif (~isnumeric(parameters.mw_pwr))||(~isreal(parameters.mw_pwr))||(numel(parameters.mw_pwr)~=1)
    error('parameters.mw_pwr should be a real number.');
end
if ~isfield(parameters,'mw_frq')
    error('microwave frequency should be specified in parameters.mw_frq variable.');
elseif (~isnumeric(parameters.mw_frq))||(~isreal(parameters.mw_frq))||(numel(parameters.mw_frq)~=1)
    error('parameters.mw_freq must be a real number.');
end
if ~isfield(parameters,'rho0')
    error('thermal equilibrium state vector must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('one or more detection states must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'mw_oper')
    error('microwave irradiation operator must be specified in parameters.mw_oper variable.');
end
if ~isfield(parameters,'ez_oper')
    error('electron Lz operator must be specified in parameters.ez_oper variable.');
end
end

% Great art is the contempt of a great man for small art.
%
% F. Scott Fitzgerald

