% Converts serpentine indexing of bosonic monomials: 
%
%                (0,0)(0,1)(0,2)     (0)(2)(5)
%                (1,0)(1,1)(1,2) <=> (1)(4)(7)
%                (2,0)(2,1)(2,2)     (3)(6)(8)
%
% into creation and annihilation operator powers:
%
%                    B(k,q)=(Cr^k)*(An^q)
%
% Syntax: 
%
%                      [K,Q]=lin2kq(N,I)
%
% Parameters:
%
%       N   - matrix dimension
%
%       I   - linear serpentine index with I=0 
%             corresponding to K=0, Q=0; an ar-
%             ray of non-negative integers 
%
% Outputs:
%
%       K   - creation operator power, an array
%             of the same dimensions as I
%
%       Q   - annihilation operator power, an ar-
%             ray of the same dimensions as I
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lin2kq.m>

function [K,Q]=lin2kq(N,I)

% Check consistency
grumble(N,I);

% Sepentine matrix
S=serpentine(N);

% Direct look-up
K=zeros(size(I));
Q=zeros(size(I));
parfor n=1:numel(I)
    [K(n),Q(n)]=find(S==I(n));
end

% Zero base
K=K-1; Q=Q-1;

end

% Consistency enforcement
function grumble(N,I)
if ~isscalar(N), error('N must be a scalar.'); end
if (~isnumeric(N))||(~isreal(N))||(mod(N,1)~=0)||...
   (~isnumeric(I))||(~isreal(I))||any(mod(I,1)~=0,'all')||any(I<0,'all')
    error('elements of I must be non-negative integers.');
end
if N<1, error('N must be a positive real integer.'); end
if any((I+1)>N^2,'all')
    error('I overflows matrix dimension.');
end
end

% "Sulphur Philosophorum": 
%      God knows what the Chymists mean by it.
%
% 1657 Physical Dictionary

