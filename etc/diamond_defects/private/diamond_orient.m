% Crystal orientation matrix for diamond defects. Syntax:
%
%          C=diamond_orient(orientation)
%
% This internal helper returns a crystal-to-laboratory rotation matrix
% for diamond-defect constructors.
%
% Parameters:
%
%    orientation  - '111', '110', or '100' crystal plane normal aligned
%                   with the laboratory z axis
%
% Outputs:
%
%    C            - crystal-to-laboratory rotation matrix
%
% <https://spindynamics.org/wiki/index.php?title=diamond_orient.m>

function C=diamond_orient(orientation)

% Check consistency
grumble(orientation);

% Select the crystal plane normal
switch orientation
    case '111'
        C=rotmat_align([1 1 1],[0 0 1]);
    case '110'
        C=rotmat_align([1 1 0],[0 0 1]);
    case '100'
        C=rotmat_align([1 0 0],[0 0 1]);
    otherwise
        error('unknown orientation specification.');
end

end

% Consistency enforcement
function grumble(orientation)
if ~ischar(orientation)
    error('orientation must be a character string.');
end
end

% A rotation is only useful when its convention is visible.

