% Converts linear indexing of bosonic monomials into K,Q indexing.
% In the linear indexing convention, bosonic monomials are listed
% in the order of increasing overall power K, and, within powers,
% in the order of increasing Q in (Cr^(Q))*(An^(K-Q)) definition.
% Zero base counting is used: 
%             
%                  I=0 -> (K=0,Q=0)
%                  I=1 -> (K=1,Q=0)
%                  I=2 -> (K=1,Q=1), et cetera...
% 
% Syntax: 
%
%                         [K,Q]=lin2kq(I)
%
% Parameters:
%
%       I   - linear indices of bosonic monomials,
%             with I=0 corresponding to K=0, Q=0;
%             an array of any dimensions
%
% Outputs:
%
%       K   - powers of bosonic monomials, an ar-
%             ray of the same dimensions as I
%
%       Q   - Q in (Cr^(Q))*(An^(K-Q)), an array
%             of the same dimensions as I
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lin2kq.m>

function [K,Q]=lin2kq(I)

% Check consistency
grumble(I);

% Get power and index
K=floor((sqrt(8*I+1)-1)/2); 
Q=I-K*(K+1)/2;

% Make sure the conversion is correct
if nnz(kq2lin(K,Q)~=I)>0
    error('IEEE arithmetic breakdown, please contact the developer.');
end

end

% Consistency enforcement
function grumble(I)
if (~isnumeric(I))||(~isreal(I))||any(mod(I,1)~=0,'all')||any(I<0,'all')
    error('all elements of the input array must be non-negative integers.');
end
end

% "Sulphur Philosophorum": 
%      God knows what the Chymists mean by it.
%
% 1657 Physical Dictionary

