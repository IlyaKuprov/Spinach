% Multiplies an RCV sparse matrix by a numeric scalar,
% in either operand order. Syntax:
%
%                        C=times(A,B)
%
% Parameters:
%
%    A,B - an RCV sparse matrix and a numeric
%          scalar, in either order
%
% Outputs:
%
%    C - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/times.m>

function C=times(A,B)

% Check consistency
grumble(A,B);

% RCV sparse by a scalar
if isa(A,'rcv')

    A.val=A.val*B; C=A;

% Scalar by RCV sparse
else

    B.val=B.val*A; C=B;

end

end

% Consistency enforcement
function grumble(A,B)
if ~xor(isa(A,'rcv'),isa(B,'rcv'))
    error('exactly one of the arguments must be an RCV sparse matrix.');
end
if isa(A,'rcv')
    scalar_arg=B;
else
    scalar_arg=A;
end
if (~isnumeric(scalar_arg))||(~isscalar(scalar_arg))
    error('multiplication is only defined by numeric scalars.');
end
end

% They say that the fish that gets away
% looks bigger than it really is.
%
% Seven Samurai film (1954)

