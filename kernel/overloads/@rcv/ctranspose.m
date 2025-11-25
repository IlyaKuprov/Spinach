% Returns the conjugate transpose of an RCV sparse matrix. Syntax:
%
%                        A=ctranspose(A)
%
% Parameters:
%
%    A   - an RCV sparse matrix
%
% Outputs:
%
%    A   - an RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/ctranspose.m>

function A=ctranspose(A)

% Check consistency
grumble(A);

% Efficiently Swap rows and columns
[A.col,A.row]=deal(A.row,A.col);

% Conjugate the values
A.val=conj(A.val);

end

% Consistency enforcement
function grumble(A)
if ~isa(A,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% At the Tower of London we remembered in
% our prayers the elephant kept there by
% James I, which, poor creature, was never
% given anything to drink but wine.
%
% Tom Holland

