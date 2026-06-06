% Converts a Hamiltonian action object into a sparse matrix. Syntax:
%
%                         S=sparse(H)
%
% Parameters:
%
%     H - a Hamiltonian action object
%
% Outputs:
%
%     S - a sparse matrix with the same action
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian_action/sparse.m>

function S=sparse(H)

% Start with a sparse matrix
S=spalloc(H.dim,H.dim,0);

% Loop over descriptor terms
for n=1:height(H.descr)

    % Compute forward and adjoint coefficients
    coeff_fwd=H.coeff_iso(n)+H.coeff_aniso(n)/2;
    coeff_adj=conj(H.coeff_aniso(n))/2;

    % Add the forward term
    if abs(coeff_fwd)>H.zero_tol
        xyz=term_oper(H,n);
        S=S+sparse(xyz(:,1),xyz(:,2),coeff_fwd*xyz(:,3),H.dim,H.dim);
    end

    % Add the adjoint anisotropic term
    if abs(coeff_adj)>H.zero_tol
        xyz=term_oper(H,n);
        S=S+sparse(xyz(:,2),xyz(:,1),coeff_adj*conj(xyz(:,3)),H.dim,H.dim);
    end

end

end

