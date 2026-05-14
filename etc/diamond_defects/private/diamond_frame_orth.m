% Right-handed frame orthogonalisation. Syntax:
%
%          frame=diamond_frame_orth(frame)
%
% This internal helper orthogonalises a supplied frame and enforces
% right-handed orientation.
%
% Parameters:
%
%    frame  - three approximate axes as columns
%
% Outputs:
%
%    frame  - orthonormal right-handed frame
%
% <https://spindynamics.org/wiki/index.php?title=diamond_frame_orth.m>

function frame=diamond_frame_orth(frame)

% Check consistency
grumble(frame);

% Orthogonalise the frame
[frame,~]=qr(frame,0);

% Enforce right-handed orientation
if det(frame)<0
    frame(:,3)=-frame(:,3);
end

end

% Consistency enforcement
function grumble(frame)
if(~isnumeric(frame)||~isreal(frame)||~isequal(size(frame),[3 3])||rank(frame)<3)
    error('frame must be a full-rank real 3x3 matrix.');
end
end

% Right-handed frames keep signs from becoming hidden parameters.

