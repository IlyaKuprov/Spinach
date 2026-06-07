% Van Vleck perturbation theory, following Shavitt and Redmon, but
% excluding the quasi-degenerate split. Syntax:
%
%                     [Ep,G]=vvpert(E0,H1,order)
%
% Parameters:
%
%     E0    - eigenvalues of H0, a column vector of real 
%             numbers
%
%     H1    - perturbation, written in the basis that di-
%             agonalises H0
%
%     order - order of perturbation theory to be used; numerical
%             artefacts appear beyond about 10-12 for typical
%             problems
%
% Outputs:
%
%     Ep    - eigenvalues of H0+H1 to the specified order,
%             a column vector of reals, not necessarily 
%             sorted in the same way as the input
%
%     G     - Van Vleck transformation generator, such that
%             expm(G) is a square unitary matrix with eigen-
%             vectors in columns, in the same order as the 
%             eigenvalues in Ep
%
% Notes: there must be no degeneracies in H0; H1 must be Hermitian,
%        the theory only converges when norm(H1,2) is much smaller
%        than the smallest energy gap in H0; complexity is cubic
%        both in the order and in the matrix dimension.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=vvpert.m>

function [Ep,G]=vvpert(E0,H1,order)

% Check consistency
grumble(E0,H1,order);

% Reciprocal energy differences
Q=1./(E0'-E0); Q(logical(speye(size(Q))))=0;
if any(~isfinite(Q),'all')
    error('H0 has degenerate energy levels.');
end

% Energy differences
delta=E0-E0.';

% Baker-Campbell-Hausdorff coefficients
bch=1./factorial(0:order);

% Allocate recursion storage
G=cell(order,1); W=cell(order,1);
C=cell(order+1,order+1);

% Compute each perturbation order
for n=1:order

    % Build higher nested commutators independent of the current generator
    for m=2:n
        K=zeros(size(H1),'like',H1);
        for k=(m-1):(n-1)
            K=K+comm(C{m,k+1},G{n-k});
        end
        C{m+1,n+1}=K;
    end

    % Build the transformed Hamiltonian term without the H0 commutator
    if n==1
        K=H1;
    else
        K=comm(H1,G{n-1});
    end
    for m=2:n
        K=K+bch(m+1)*C{m+1,n+1};
    end

    % Split the term into the effective Hamiltonian and generator
    W{n}=diag(diag(K));
    G{n}=Q.*(K-W{n});

    % Store the first nested commutator for future orders
    C{2,n+1}=delta.*G{n};
    if n>1
        C{2,n+1}=C{2,n+1}+comm(H1,G{n-1});
    end

end

% Sum the energy corrections
Ep=E0;
for n=1:order
    Ep=Ep+real(diag(W{n}));
end

% Sum the generator corrections
G=reshape(G,[1 1 numel(G)]);
G=sum(cell2mat(G),3);

end

% Consistency enforcement
function grumble(E0,H1,order)
if (~isnumeric(E0))||(~isreal(E0))||(~iscolumn(E0))
    error('E0 must be a real column vector.');
end
if (~isnumeric(H1))||(~ishermitian(H1))
    error('H1 must be a Hermitian matrix.');
end
if (numel(E0)~=size(H1,1))||(numel(E0)~=size(H1,2))
    error('dimensions of E0 and H1 must be consistent.');
end
if (~isnumeric(order))||(~isreal(order))||(~isscalar(order))||...
   (mod(order,1)~=0)||(order<1)
    error('order must be a positive real integer.');
end
end

% When have we ever had stability and security in this
% region that we should be concerned about losing it?
%
% Kurdish officials, to BBC

