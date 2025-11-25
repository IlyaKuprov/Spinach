% Converts an RCV sparse matrix into a Matlab sparse 
% matrix. Syntax:
%                      A=sparse(A)
%
% Parameters:
%
%    A   - RCV sparse matrix
%
% Outputs:
%
%    A   - Matlab sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/sparse.m>

function A=sparse(A)

% Check consistency
grumble(A);

% Check if empty
if isempty(A.col)

    % Empty matrix of a specified size
    A=spalloc(A.numRows,A.numCols,0);

else

    % Call Matlab's sparse matrix constructor
    A=sparse(A.row,A.col,A.val,A.numRows,A.numCols);

end

end

% Consistency enforcement
function grumble(A)
if ~isa(A,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% Working 16 hours a day, 7 days a week, 52 weeks
% in a year, and people still calling me lucky.
%
% Elon Musk

