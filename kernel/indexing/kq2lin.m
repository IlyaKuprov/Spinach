% Converts k,q indexing of single transition operators into 
% serpentine indexing:
%
%               (1,1)(1,2)(1,3)     (1)(3)(6)
%               (2,1)(2,2)(2,3) <=> (2)(5)(8)
%               (3,1)(3,2)(3,3)     (4)(7)(9)
%
% Syntax: 
%
%                       I=kq2lin(N,K,Q)
%
% Parameters:
%
%       N   - matrix dimension, a scalar
%
%       K   - first index, positive integer
%             array of any size
%
%       Q   - second index, positive integer
%             array of the same size as K
%
% Outputs:
%
%       I   - linear serpentine index, an ar-
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
    I(n)=S(K(n),Q(n));
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
if any(K<1,'all')||any(Q<1,'all')
    error('elements of K and Q must be positive.');
end
if any(size(K)~=size(Q),'all')
    error('K and Q arrays must have the same size.');
end
if any(K>N,'all')||any(Q>N,'all')
    error('K,Q overflow matrix dimension.');
end
end

% Crisis. Survival. Advancement.
%
% Frank Herbert, in "Dune: Prophesy"

