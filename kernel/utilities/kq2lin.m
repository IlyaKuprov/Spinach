% Converts K,Q indexing of bosonic monomials into linear indexing.
% In the linear indexing convention, bosonic monomials are listed
% in the order of increasing overall power K, and, within powers,
% in the order of increasing Q in (Cr^(Q))*(An^(K-Q)) definition.
% Zero base counting is used: 
%
%                  (K=0,Q=0) -> I=0
%                  (K=1,Q=0) -> I=1
%                  (K=1,Q=0) -> I=2, et cetera...
%
% Syntax: 
%
%                         I=kq2lin(K,Q)
%
% Parameters:
%
%       K   - powers of bosonic monomials, an ar-
%             ray of any dimensions
%
%       Q   - Q in (Cr^(Q))*(An^(K-Q)), an array
%             of the same dimensions as K
%
% Outputs:
%
%       I   - linear indices of bosonic monomials,
%             with I=0 corresponding to K=0, Q=0;
%             array of the same dimensions as K
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kq2lin.m>

function I=kq2lin(K,Q)

% Check consistency
grumble(K,Q);

% Get the linear index
I=K.*(K+1)/2+Q;

end

% Consistency enforcement
function grumble(K,Q)
if (~isnumeric(K))||(~isreal(K))||any(mod(K,1)~=0,'all')||...
   (~isnumeric(Q))||(~isreal(Q))||any(mod(Q,1)~=0,'all')
    error('all elements of the inputs must be real integers.');
end
if any(K<0,'all')
    error('bosonic monomial power must be non-negative.');
end
if any(size(K)~=size(Q),'all')
    error('array dimensions are inconsistent.');
end
if any(Q<0,'all')||any(Q>K,'all')
    error('Q index must be between 0 and K.');
end
end

% Crisis. Survival. Advancement.
%
% Frank Herbert, in "Dune: Prophesy"

