% Matrix products involving a Hamiltonian action object. Syntax:
%
%                         c=mtimes(a,b)
%
% Parameters:
%
%     a,b   - Hamiltonian action object and a numeric array
%
% Outputs:
%
%     c     - multiplication result
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian_action/mtimes.m>

function c=mtimes(a,b)

% Process Hamiltonian action on a numeric array
if isa(a,'hamiltonian_action')&&isnumeric(b)

    % Check matrix dimensions
    if size(b,1)~=size(a,2)
        error('matrix dimensions are not consistent.');
    end

    % Preallocate the answer
    c=complex(zeros(a.dim,size(b,2),'like',b));

    % Loop over descriptor terms
    for n=1:height(a.descr)

        % Compute forward and adjoint coefficients
        coeff_fwd=a.coeff_iso(n)+a.coeff_aniso(n)/2;
        coeff_adj=conj(a.coeff_aniso(n))/2;

        % Apply the forward term
        if abs(coeff_fwd)>a.zero_tol
            xyz=term_oper(a,n);
            c=c+sparse(xyz(:,1),xyz(:,2),coeff_fwd*xyz(:,3),a.dim,a.dim)*b;
        end

        % Apply the adjoint anisotropic term
        if abs(coeff_adj)>a.zero_tol
            xyz=term_oper(a,n);
            c=c+sparse(xyz(:,2),xyz(:,1),coeff_adj*conj(xyz(:,3)),a.dim,a.dim)*b;
        end

    end

    % Add the oriented giant spin contribution
    if nnz(a.giant)>0
        c=c+a.giant*b;
    end

    % Return the result
    return

end

% Complain and bomb out
error('left operand must be a Hamiltonian action object and right operand must be numeric.');

end

