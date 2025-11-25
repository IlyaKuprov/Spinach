% Horizontal concatenation for RCV sparse matrices. Syntax:
%
%                        A=horzcat(A,B)
%
% Parameters:
%
%    A   - left RCV sparse matrix
%
%    B   - right RCV sparse matrix
%
% Outputs:
%
%    A   - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/horzcat.m>

function A=horzcat(A,B)

% Check consistency
grumble(A,B);

% Align locations
if A.isGPU||B.isGPU
    A=gpuArray(A);
    B=gpuArray(B);
end

% Shift column indices
B.col=B.col+A.numCols;

% Concatenate RCV arrays
A.row=[A.row; B.row];
A.col=[A.col; B.col];
A.val=[A.val; B.val];

% Update column count
A.numCols=A.numCols+B.numCols;

end

% Consistency enforcement
function grumble(A,B)
if ~isa(A,'rcv')||~isa(B,'rcv')
    error('both inputs must be RCV sparse matrices.');
end
if A.numRows~=B.numRows
    error('row counts must match for horizontal concatenation.');
end
end

% The back half of your forties is a cursed age. It's not
% so much that (as the cliché goes) the policemen look
% younger - it's more that you're gloomily aware that
% you haven't done anything that they'd even consider
% arresting you for.
%
% Sam Leith

