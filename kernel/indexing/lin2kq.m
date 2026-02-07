% Converts serpentine indices of single transition operators: 
%
%                (1,1)(1,2)(1,3)     (1)(3)(6)
%                (2,1)(2,2)(2,3) <=> (2)(5)(8)
%                (3,1)(3,2)(3,3)     (4)(7)(9)
%
% into their non-zero element indices. Syntax: 
%
%                      [K,Q]=lin2kq(N,I)
%
% Parameters:
%
%       N   - matrix dimension
%
%       I   - linear serpentine index, an ar-
%             ray of positive integers 
%
% Outputs:
%
%       K   - row index of the non-zero element,
%             an array of the same size as I
%
%       Q   - col index of the non-zero element,
%             an array of the same size as I
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

end

% Consistency enforcement
function grumble(N,I)
if ~isscalar(N), error('N must be a scalar.'); end
if (~isnumeric(N))||(~isreal(N))||(mod(N,1)~=0)||...
   (~isnumeric(I))||(~isreal(I))||any(mod(I,1)~=0,'all')||any(I<1,'all')
    error('elements of I must be positive integers.');
end
if N<1, error('N must be a positive real integer.'); end
if any(I>N^2,'all')
    error('I overflows matrix dimension.');
end
end

% "Sulphur Philosophorum": 
%      God knows what the Chymists mean by it.
%
% 1657 Physical Dictionary

