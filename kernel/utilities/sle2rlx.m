% Relaxation superoperator extracted from a stochastic Liouville equation
% generator by Markovian back-projection onto the spin Liouville space.
% The SLE generator acts on the direct product of a lattice space (Wigner
% functions of molecular orientation in gridfree.m, or chemical species)
% and the spin Liouville space; the lattice dynamics is projected out and
% appears in the spin space as a Redfield-type relaxation superoperator.
% Syntax:
%
%                 [R,L0]=sle2rlx(spin_system,L,pi_vec,one_vec)
%
% Parameters:
%
%    L        - SLE generator, L=H+1i*R+1i*K, in the kron(lattice,spin)
%               index order used by gridfree.m, sparse or polyadic
%
%    pi_vec   - lattice equilibrium distribution (right null vector of
%               the lattice generator), column vector
%
%    one_vec  - lattice integration functional (left null vector of
%               the lattice generator), column vector such that the
%               contraction is one_vec'*(...); one_vec'*pi_vec must
%               be non-zero, the normalisation is done internally
%
% Outputs:
%
%    R        - relaxation superoperator in the spin Liouville space,
%               lab frame, in the convention of relaxation.m (it en-
%               ters the Liouvillian as 1i*R, dynamic frequency shi-
%               fts are its imaginary part)
%
%    L0       - lattice-averaged Liouvillian <one|L|pi>
%
% Note: with P=|pi><one|(x)1 and Q=1-P, Nakajima-Zwanzig partitioning
%       gives the exact spin-space resolvent [s-L0-C*(s-A)^(-1)*B]^(-1)
%       with A=QLQ, B=QL|pi>, and C=<one|LQ. The Markov limit is the
%       Sylvester equation A*X-X*L0=-B, and then 1i*R=C*X. The columns
%       of X are obtained from the Schur form of L0 by sparse solves
%       of bordered systems that keep the lattice equilibrium compo-
%       nent of X at zero, so pi_vec and one_vec need not be sparse.
%
% Note: the result is a valid relaxation superoperator when the lattice
%       relaxes faster than the spins it modulates (Redfield condition);
%       the eigenvalues of L0+1i*R should then approximate the spin-like
%       branch of the SLE spectrum. The function does not check this.
%
% talos@spindynamics.org
%
% <https://spindynamics.org/wiki/index.php?title=sle2rlx.m>

function [R,L0]=sle2rlx(spin_system,L,pi_vec,one_vec)

% Check consistency
grumble(L,pi_vec,one_vec);

% Inflate polyadic generators
L=inflate(L);

% Lattice and spin space dimensions
spc_dim=numel(pi_vec); spn_dim=size(L,1)/spc_dim;
report(spin_system,['lattice space dimension ' num2str(spc_dim) ...
                    ', spin space dimension ' num2str(spn_dim)]);

% Normalise the lattice equilibrium
pi_vec=pi_vec/(one_vec'*pi_vec);

% Injection into and contraction out of the lattice equilibrium
P_inj=kron(sparse(pi_vec),speye(spn_dim));
P_con=kron(sparse(one_vec'),speye(spn_dim));

% Lattice-averaged Liouvillian
L0=P_con*L*P_inj;

% Schur form of the lattice-averaged Liouvillian
[U,T]=schur(full(L0),'complex');

% Solve the Sylvester equation column by column
report(spin_system,'solving the Sylvester equation...');
X=zeros(spc_dim*spn_dim,spn_dim); Z=sparse(spn_dim,spn_dim);
for n=1:spn_dim
    rhs=-L*(P_inj*U(:,n))+X(:,1:(n-1))*T(1:(n-1),n);
    sol=[L-T(n,n)*speye(spc_dim*spn_dim), -P_inj;
         P_con,                           Z]\[rhs; zeros(spn_dim,1)];
    X(:,n)=sol(1:(spc_dim*spn_dim));
end
X=X*U';

% Relaxation superoperator
R=-1i*P_con*L*X;

% Report the outcome
report(spin_system,['relaxation superoperator norm ' num2str(norm(R,1))]);

end

% Consistency enforcement
function grumble(L,pi_vec,one_vec)
if (~isnumeric(L))&&(~isa(L,'polyadic'))
    error('L must be a numeric or polyadic array.');
end
if size(L,1)~=size(L,2)
    error('L must be a square matrix.');
end
if (~isnumeric(pi_vec))||(~iscolumn(pi_vec))||...
   (~isnumeric(one_vec))||(~iscolumn(one_vec))
    error('pi_vec and one_vec must be numeric column vectors.');
end
if numel(pi_vec)~=numel(one_vec)
    error('pi_vec and one_vec must have the same length.');
end
if mod(size(L,1),numel(pi_vec))~=0
    error('the dimension of L must be a multiple of the lattice dimension.');
end
if one_vec'*pi_vec==0
    error('one_vec''*pi_vec must be non-zero.');
end
end

% There is no such thing as a free lunch, but there is such a
% thing as a lunch that somebody else has already paid for.
%
% Anonymous

