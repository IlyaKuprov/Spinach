% Cross-polarisation experiment in the rotating frame. Applies an
% ideal pi/2 pulse using the specified operators, then evolves the
% system with the specified spin-lock terms added to the Liovilli-
% an. The contact curve is returned. Syntax:
%
%   contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.irr_powers - a matrix containing the values
%                             of the spin-lock nutation fre-
%                             quency on each channel (rows)
%                             at each time slice (cols), Hz
%                            
%     parameters.irr_opers  - a cell array of spin operators
%                             corresponding to the spin-lock
%                             on each channel
%
%     parameters.exc_opers  - a cell array of spin operators
%                             for the ideal pi/2 excitation 
%                             pulse (same flip angle on all
%                             channels)
%
%     parameters.time_steps - a vector of time slice durati-
%                             ons, seconds
%
%     parameters.rho0       - initial state vector
%
%     parameters.coil       - detection state vector
%
%     H - Hamiltonian matrix, received from context function
%
%     R - relaxation superoperator, received from context function
%
%     K - kinetics superoperator, received from context function
%
% Output:
%
%     contact_curve - contact curve detected on the coil
%                     state specified in parameters.coil
%
% marina.carravetta@soton.ac.uk
% p.t.williamson@soton.ac.uk
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cp_contact_hard.m>

function contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Build and project the excitation operator
exc_oper=parameters.exc_opers{1};
for n=2:numel(parameters.exc_opers)
    exc_oper=exc_oper+parameters.exc_opers{n};
end
exc_oper=kron(speye(parameters.spc_dim),exc_oper);

% Apply a perfect hard pi/2 excitation pulse
rho=step(spin_system,exc_oper,parameters.rho0,pi/2);

% Preallocate contact curves
contact_curve=zeros(size(parameters.coil,2),...
                    numel(parameters.time_steps)+1);
contact_curve(:,1)=parameters.coil'*rho;

% Loop over the time steps
for n=1:numel(parameters.time_steps)

    % Build and project the spin-lock operator
    irr_oper=2*pi*parameters.irr_powers(1,n)*...
                  parameters.irr_opers{1};
    for k=2:numel(parameters.irr_opers)
        irr_oper=irr_oper+2*pi*parameters.irr_powers(k,n)*...
                               parameters.irr_opers{k};
    end
    irr_oper=kron(speye(parameters.spc_dim),irr_oper);
 
    % Apply the evolution step
    rho=step(spin_system,L+irr_oper,rho,parameters.time_steps(n));

    % Project out the observable
    contact_curve(:,n+1)=parameters.coil'*rho;

end

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'irr_powers')
    error('nutation frequencies must be specified in parameters.irr_powers variable.');
end
if ~isfield(parameters,'irr_opers')
    error('spin-lock operators must be specified in parameters.irr_opers variable.');
end
if ~iscell(parameters.irr_opers)
    error('parameters.irr_opers must be a cell array of matrices.');
end
if ~isfield(parameters,'exc_opers')
    error('excitation operators must be specified in parameters.exc_opers variable.');
end
if ~iscell(parameters.exc_opers)
    error('parameters.exc_opers must be a cell array of matrices.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
end

% No one would remember the Good Samaritan if he'd
% only had good intentions; he had money as well.
%
% Margaret Thatcher

