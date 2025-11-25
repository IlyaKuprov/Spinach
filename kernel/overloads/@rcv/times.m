% Multiplies an RCV sparse matrix by a numeric scalar. Syntax:
%
%                         A=times(k,A)
%
% Parameters:
%
%    k - numeric scalar
%
%    A - RCV sparse matrix
%
% Outputs:
%
%    A - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/times.m>

function A=times(k,A)

% Check consistency
grumble(k,A);

% Multiply values
A.val=k*A.val;

end

% Consistency enforcement
function grumble(k,A)
if ~isa(A,'rcv')
    error('the second argument must be an RCV sparse matrix.');
end
if (~isnumeric(k))||(~isscalar(k))
    error('multiplication is only defined by numeric scalars.');
end
end

% They say that the fish that gets away
% looks bigger than it really is.
%
% Seven Samurai film (1954)

