% Eigensystem of sparse Hamiltonians to user-specified order in
% RSPT with careful handling of diagonal dominanace and an opti-
% on to do exact diagonalisation (expensive). Syntax:
%
%                     [E,V]=rspt_eig(H,order)
%
% Parameters:
%
%    H     - sparse Hermitian matrix; non-Hermitian parts 
%            will be dropped, perturbation theory will be
%            done with respect to the off-diagonal part
%
%    order - RSPT order, use Inf for exact diagonalisation
%
% Outputs:
%
%    E     - a column vector of energies, sorted in ascen-
%            ding order
%
%    V     - a matrix with eigenvectors in columns, sorted 
%            left to right in the same order as the energies
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rspt_eig.m>

function [E,V]=rspt_eig(H,order)

% Decide the method
switch order

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

% Eigensystem sorting
[E,idx]=sort(E,'ascend'); V=V(:,idx);

end

% Мне бы водки речушку
% Да баб деревеньку.
% Я бы пил потихоньку
% И е..ал помаленьку.
% 
% Сергей Есенин

