% Truncates the singular value decomposition at the specified rank
% and reassembles the matrix. Syntax:
%
%                        A=keep_rank(A,rank)
%
% Parameters:
%
%         A     -  real or complex matrix, will be 
%                  converted to full if a sparse
%                  matrix is received
%
%         nsvk  -  number of singular values to keep
%
% Outputs:
%
%         A     -  filtered matrix, returned as full
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=keep_rank.m>

function A=keep_rank(A,nsvk)

% Check consistency
grumble(A,nsvk);

% Run singular value decomposition
[U,S,V]=svd(full(A));

% Truncate to the specified rank and rebuild
A=U(:,1:nsvk)*S(1:nsvk,1:nsvk)*V(:,1:nsvk)';

end

% Consistency enforcement
function grumble(A,nsvk)
if (~isnumeric(A))||(size(A,1)<=1)||(size(A,2)<=1)
    error('A must be a matrix.');
end
if (~isnumeric(nsvk))||(~isreal(nsvk))||...
   (~isscalar(nsvk))||(nsvk<1)||...
   (mod(nsvk,1)~=0)||any(nsvk>size(A))
    error('nsvk must be a positive integer smaller than dim(A)');
end
end

% Endure. In enduring, grow strong.
% 
% Dakkon

