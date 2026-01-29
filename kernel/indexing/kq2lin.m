% Converts k,q power indexing of bosonic monomials: 
%
%                    B(k,q)=(Cr^k)*(An^q)
%
% into serpentine indexing:
%
%                (0,0)(0,1)(0,2)     (0)(2)(5)
%                (1,0)(1,1)(1,2) <=> (1)(4)(7)
%                (2,0)(2,1)(2,2)     (3)(6)(8)
%
% Syntax: 
%
%                       I=kq2lin(N,K,Q)
%
% Parameters:
%
%       N   - matrix dimension, a scalar
%
%       K   - creation operator power, an
%             integer array of any size
%
%       Q   - annihilation operator power,
%             an integer array of the same
%             size as K
%
% Outputs:
%
%       I   - linear serpentine index with I=0 
%             corresponding to K=0, Q=0; an ar-
%             ray of the same size as inputs
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kq2lin.m>

function I=kq2lin(N,K,Q)

% Check consistency
grumble(N,K,Q);

% Sepentine matrix
S=serpentine(N);

% Direct look-up
I=zeros(size(K));
for n=1:numel(K)
    I(n)=S(K(n)+1,Q(n)+1);
end

end

% Consistency enforcement
function grumble(N,K,Q)
if ~isscalar(N), error('N must be a scalar.'); end
if (~isnumeric(N))||(~isreal(N))||(mod(N,1)~=0)||...
   (~isnumeric(K))||(~isreal(K))||any(mod(K,1)~=0,'all')||...
   (~isnumeric(Q))||(~isreal(Q))||any(mod(Q,1)~=0,'all')
    error('all elements of the inputs must be real integers.');
end
if N<1, error('N must be a positive real integer.'); end
if any(K<0,'all')||any(Q<0,'all')
    error('elements of K and Q must be non-negative.');
end
if any(size(K)~=size(Q),'all')
    error('K and Q arrays must have the same size.');
end
if any((K+1)>N,'all')||any((Q+1)>N,'all')
    error('K,Q overflow matrix dimension.');
end
end

% Crisis. Survival. Advancement.
%
% Frank Herbert, in "Dune: Prophesy"

