% Eigensystem of sparse Hamiltonians to user-specified order in
% RSPT with careful handling of diagonal dominanace and an opti-
% on to do exact diagonalisation (expensive). Syntax:
%
% [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,alpha)
%
% Parameters:
%
%     Hz    -  laboratory frame Hamiltonian, containing only Zee-
%              man terms at 1 Tesla
%
%     Hc    -  laboratory frame Hamiltonian, containing all spin-
%              spin couplings, but no Zeeman terms
%
%     Hmw   -  observable operator (Hilbert space) or observable
%              vector (Liouville space), without the amplitude 
%              prefactor
%
%     alpha - magnetic field, Tesla
%
%     parameters.rspt_order - perturbation theory order to use
%                             to account for the off-diagonal
%                             part of the Hamiltonian, Inf for
%                             exact diagonalisation
%
% Outputs:
%
%    E     - a column vector of energies, sorted in ascen-
%            ding order
%
%    V     - a matrix with eigenvectors in columns, sorted 
%            left to right in the same order as the energies
%
%    dE    - a column vector of dE/dB derivatives, sorted in
%            the same order as the energies
%
%    T     - a matrix of transition moments under Hmw
%
%    LP    - a column vector of energy level populations, sor-
%            ted in the same order as the energies
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rspt_eig.m>

function [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,alpha)

% Recursive call for symmetry
if isfield(spin_system.bas,'irrep')

    % Preallocate irrep output blocks
    n_irreps=numel(spin_system.bas.irrep);
    E=cell(n_irreps,1);  V=cell(1,n_irreps); dE=cell(n_irreps,1);

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
        [E{n},V{n},dE{n}]=rspt_eig(spin_system_nosym,parameters,...
                                   HzIrr,HcIrr,HmwIrr,alpha);

        % Project back
        V{n}=P*V{n};

    end

    % Concatenate irrep blocks
    E=cell2mat(E); V=cell2mat(V); dE=cell2mat(dE);

    % Sort the energies in ascending order
    [E,idx]=sort(E,'ascend'); V=V(:,idx); dE=dE(idx);

else

    % Hamiltonian
    H=alpha*Hz+Hc;

    % Decide the method
    switch parameters.rspt_order

        case {1,2,3,4}

            % Diag and off-diag Hamiltonian
            H_diag=diag(H); H_offd=(H+H')/2;
            H_offd(logical(speye(size(H))))=0;

            % Call perturbation theory
            [E,V]=rspert(H_diag,H_offd,order);

        case Inf

            % Full diagonalisation
            [V,E]=eig(full((H+H')/2),'vector');

        otherwise

            % Complain and bomb out
            error('unsupported perturbation theory order.');

    end

    % Sort the energies in ascending order
    [E,idx]=sort(E,'ascend'); V=V(:,idx);

    % Hellmann-Feynman derivatives, if needed
    if nargout>2, dE=real(diag(V'*Hz*V)); end

end

% Transition moments, if required
if nargout>3, T=abs(V'*Hmw*V).^2; end

% Energy level populations, if required
if (nargout>4)&&isfield(parameters,'rho0')

    % User-specified density matrix
    LP=real(diag(V'*parameters.rho0*V));

elseif nargout>4

    % Thermal equilibrium
    rho=equilibrium(spin_system,H);
    LP=real(diag(V'*rho*V));

end

end

% Мне бы водки речушку
% Да баб деревеньку.
% Я бы пил потихоньку
% И е..ал помаленьку.
% 
% Сергей Есенин

