% Converts an RCV sparse matrix into a full matrix. Syntax:
%
%                       A=full(A)
%
% Parameters:
%
%    A   - an RCV sparse matrix
%
% Outputs:
%
%    A   - a full Matlab matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/full.m>

function A=full(A)

% Check consistency
grumble(A);

% Delegate to Matlab
A=full(sparse(A));

end

% Consistency enforcement
function grumble(A)
if ~isa(A,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% Whenever you find yourself on the side of the
% majority, it is time to pause and reflect.
%
% Mark Twain

