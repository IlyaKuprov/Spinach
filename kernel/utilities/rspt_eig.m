% Eigensystem of sparse Hamiltonians to user-specified order in
% RSPT with careful handling of diagonal dominanace and an opti-
% on to do exact diagonalisation (expensive). Syntax:
%
%  [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)
%
% Parameters:
%
%     Hz   -  laboratory frame Hamiltonian, containing only 
%             Zeeman terms at 1 Tesla
%
%     Hc   -  laboratory frame Hamiltonian, containing all 
%             spin-spin couplings, but no Zeeman terms
%
%     Hmw  -  observable operator without the amplitude pre-
%             factor (2-norm should be around 1)
%
%     B    -  magnetic field, Tesla
%
%     parameters.rspt_order - perturbation theory order to use
%                             to account for the off-diagonal
%                             part of the Hamiltonian, Inf for
%                             exact diagonalisation
%
%     parameters.rho0 - [optional] when a matrix, sets a user-
%                       specified thermal equilibrium state;
%                       when a function handle f(B,alp,bet,gam)
%                       sets the function to call to obtain
%                       the thermal equilibrium at each orien-
%                       tation and magnetic field; if not pro-
%                       vided, the thermal equilibrium is com-
%                       puted at the current temperature, ori-
%                       entation, and magnetic field
%
% Outputs:
%
%     E    - a column vector of energies, sorted in ascen-
%            ding order (rad/s)
%
%     V    - a matrix with eigenvectors in columns, sorted 
%            left to right in the same order as the energies
%
%     dE   - a column vector of dE/dB derivatives, sorted in
%            the same order as the energies
%
%     T    - a matrix of transition moments under Hmw
%
%     LP   - a column vector of energy level populations, sor-
%            ted in the same order as the energies
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rspt_eig.m>

function [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)

% Check consistency
grumble(parameters,Hz,Hc,Hmw,B);

% Recursive call for symmetry
if isfield(spin_system.bas,'irrep')

    % Preallocate irrep output blocks
    n_irreps=numel(spin_system.bas.irrep);
    E=cell(n_irreps,1); V=cell(1,n_irreps);

    % Loop over irreps
    for n=1:n_irreps

        % Extract irrep projector
        P=spin_system.bas.irrep(n).projector;

        % Project Hamiltonian components
        HzIrr=P'*Hz*P; HzIrr=(HzIrr+HzIrr')/2; 
        HcIrr=P'*Hc*P; HcIrr=(HcIrr+HcIrr')/2; 

        % Project microwave operator
        HmwIrr=P'*Hmw*P; HmwIrr=(HmwIrr+HmwIrr')/2;

        % Issue a recursive call
        spin_system_nosym=spin_system;
        spin_system_nosym.bas=rmfield(spin_system.bas,'irrep');
        [E{n},V{n}]=rspt_eig(spin_system_nosym,parameters,...
                             HzIrr,HcIrr,HmwIrr,B);

        % Project back
        V{n}=P*V{n};

    end

    % Concatenate irrep blocks
    E=cell2mat(E); V=cell2mat(V);

else

    % Hamiltonian
    H=B*Hz+Hc; H=full((H+H')/2);

    % Decide the method
    switch parameters.rspt_order

        case {1,2,3,4}

            % Diag and off-diag Hamiltonian
            H_diag=diag(H); H_offd=H-diag(diag(H));
            
            % Call perturbation theory
            [E,V]=rspert(H_diag,H_offd,parameters.rspt_order);

        case Inf

            % Full diagonalisation
            [V,E]=eig(H,'vector');

        otherwise

            % Complain and bomb out
            error('unsupported perturbation theory order.');

    end

end

% Sort the energies in ascending order
[E,idx]=sort(E,'ascend'); V=V(:,idx);

% Hellmann-Feynman dE/dB, if required
if nargout>2, dE=real(diag(V'*Hz*V)); end

% Transition moments, if required
if nargout>3, T=abs(V'*Hmw*V).^2; end

% Energy level populations, if required
if (nargout>4)&&isfield(parameters,'rho0')

    % User-specified equilibrium
    if isa(parameters.rho0,'function_handle')

        % Orientation-dependent
        rho0=parameters.rho0(B,parameters.orientation(1),...
                               parameters.orientation(2),...
                               parameters.orientation(3));

    else

        % Orientation-independent
        rho0=parameters.rho0;

    end

    % Level populations
    LP=real(diag(V'*rho0*V));

elseif nargout>4

    % Thermal equilibrium
    rho0=equilibrium(spin_system,H);

    % Level populations
    LP=real(diag(V'*rho0*V));

end

end

% Consistency enforcement
function grumble(parameters,Hz,Hc,Hmw,B)
if ~isfield(parameters,'rspt_order')
    error('parameters.rspt_order subfield must be present.');
end
if (~isnumeric(parameters.rspt_order))||...
   (~isreal(parameters.rspt_order))||...
   (~isscalar(parameters.rspt_order))||...
   (parameters.rspt_order<1)
    error('parameters.rspt_order must be a positive integer or Inf.');
end
if (~isnumeric(Hz))||(size(Hz,1)~=size(Hz,2))||...
   (~isnumeric(Hc))||(size(Hc,1)~=size(Hc,2))||...
   (~isnumeric(Hmw))||(size(Hmw,1)~=size(Hmw,2))
    error('Hc, Hz, Hmw must be square matrices.');
end
if (~isnumeric(B))||(~isreal(B))||(~isscalar(B))
    error('B must be a real scalar.');
end
end

% Мне бы водки речушку
% Да баб деревеньку.
% Я бы пил потихоньку
% И е..ал помаленьку.
% 
% Сергей Есенин

