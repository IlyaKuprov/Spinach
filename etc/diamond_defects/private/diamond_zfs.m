% Zero-field tensor builder for diamond defects. Syntax:
%
%          M=diamond_zfs(D,E,frame)
%
% This internal helper builds a symmetric traceless zero-field splitting
% tensor from Spinach D and E parameters.
%
% Parameters:
%
%    D      - axial zero-field splitting parameter, Hz
%
%    E      - rhombic zero-field splitting parameter, Hz
%
%    frame  - three principal axes as columns
%
% Outputs:
%
%    M      - symmetric traceless zero-field tensor in Hz
%
% <https://spindynamics.org/wiki/index.php?title=diamond_zfs.m>

function M=diamond_zfs(D,E,frame)

% Check consistency
grumble(D,E,frame);

% Orthogonalise the principal-axis frame
frame=diamond_frame_orth(frame);

% Build the zero-field tensor
M=diamond_traceless(frame*zfs2mat(D,E,0,0,0)*frame');

end

% Consistency enforcement
function grumble(D,E,frame)
if(~isnumeric(D)||~isreal(D)||~isscalar(D))
    error('D must be a real scalar.');
end
if(~isnumeric(E)||~isreal(E)||~isscalar(E))
    error('E must be a real scalar.');
end
if(~isnumeric(frame)||~isreal(frame)||~isequal(size(frame),[3 3]))
    error('frame must be a real 3x3 matrix.');
end
end

% Quadratic splittings should stay traceless.

