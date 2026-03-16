% A simple shorthand for the anticommutator of two 
% matrices. Syntax:
%
%                   C=acomm(A,B)
%
% Parameters:
%
%   A,B - square matrices
%
% Outputs:
%
%   C   - a square matrix
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=acomm.m>

function C=acomm(A,B)

% Check consistency
grumble(A,B);

% Do the deed
C=A*B+B*A;

end

% Consistency enforcement
function grumble(A,B)
if (~isnumeric(A))||(~ismatrix(A))||(size(A,1)~=size(A,2))
    error('A must be a numeric square matrix.');
end
if (~isnumeric(B))||(~ismatrix(B))||(size(B,1)~=size(B,2))
    error('B must be a numeric square matrix.');
end
if ~isequal(size(A),size(B))
    error('A and B must have the same dimensions.');
end
end

% Enough of all this academic chatter, back
% again to devilry!
%
% Mephisto


