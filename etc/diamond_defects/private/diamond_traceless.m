% Symmetric traceless projection for diamond tensors. Syntax:
%
%          M=diamond_traceless(M)
%
% This internal helper enforces symmetric traceless form for quadratic
% spin tensors.
%
% Parameters:
%
%    M  - tensor matrix
%
% Outputs:
%
%    M  - symmetric traceless tensor matrix
%
% <https://spindynamics.org/wiki/index.php?title=diamond_traceless.m>

function M=diamond_traceless(M)

% Check consistency
grumble(M);

% Symmetrise the input tensor
M=(M+M')/2;

% Remove the isotropic part
M=M-eye(3)*trace(M)/3;

% Symmetrise numerical round-off
M=(M+M')/2;

end

% Consistency enforcement
function grumble(M)
if(~isnumeric(M)||~isreal(M)||~isequal(size(M),[3 3]))
    error('M must be a real 3x3 matrix.');
end
end

% The trace belongs in the origin, not in the splitting.

