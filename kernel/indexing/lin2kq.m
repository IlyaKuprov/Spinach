% Converts linear serpentine indexing of matrices into their
% k,q indexing. In base-1 indexing convention:
%
%               (1,1)(1,2)(1,3)     (1)(3)(6)
%               (2,1)(2,2)(2,3) <=> (2)(5)(8)
%               (3,1)(3,2)(3,3)     (4)(7)(9)
%
% and in base 0 indexing convention:
%
%               (0,0)(0,1)(0,2)     (0)(2)(5)
%               (1,0)(1,1)(1,2) <=> (1)(4)(7)
%               (2,0)(2,1)(2,2)     (3)(6)(8)
%
% Syntax: 
%
%                [K,Q]=lin2kq(N,I,idx_base)
%
% Parameters:
%
%         N   - matrix dimension
%
%         I   - linear serpentine indices, an
%               array of integers 
%
%    idx_base - indexing base, 0 or 1 
%
% Outputs:
%
%         K   - row indices, an array of the 
%               same size as I
%
%         Q   - col indices, an array of the 
%               same size as I
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lin2kq.m>

function [K,Q]=lin2kq(N,I,idx_base)

% Check consistency
grumble(N,I,idx_base);

% Sepentine matrix
S=serpentine(N,idx_base);

% Direct look-up
K=zeros(size(I));
Q=zeros(size(I));
switch idx_base

    case 0

        % 0-base indexing
        for n=1:numel(I)
            [K(n),Q(n)]=find(S==I(n));
        end
        K=K-1; Q=Q-1;

    case 1
        
        % 1-base indexing
        for n=1:numel(I)
            [K(n),Q(n)]=find(S==I(n));
        end

    otherwise

        % Complain and bomb out
        error('unsupported indexing base.');

end

end

% Consistency enforcement
function grumble(N,I,idx_base)
if (~isnumeric(idx_base))||(~isreal(idx_base))||...
   (~isscalar(idx_base))||(~ismember(idx_base,[0 1]))
    error('idx_base must be 0 or 1.');
end
if (~isnumeric(N))||(~isreal(N))||...
   (~isscalar(N))||(mod(N,1)~=0)||(N<1)
    error('N must be a positive real integer.');
end
if (~isnumeric(I))||(~isreal(I))||any(mod(I,1)~=0,'all')
    error('elements of I must be integers.');
end
if any(I<idx_base,'all')
    error('I cannot be smaller than the indexing base.');
end
if any(I>(N^2-1+idx_base),'all')
    error('element(s) of I overflow matrix dimension.');
end
end

% "Sulphur Philosophorum": 
%      God knows what the Chymists mean by it.
%
% 1657 Physical Dictionary

