% Principal-axis frame from three axes. Syntax:
%
%          frame=diamond_frame_xyz(xaxis,yaxis,zaxis)
%
% This internal helper orthogonalises three supplied axes into a
% right-handed principal-axis frame.
%
% Parameters:
%
%    xaxis  - approximate x axis
%
%    yaxis  - approximate y axis
%
%    zaxis  - approximate z axis
%
% Outputs:
%
%    frame  - orthonormal right-handed frame
%
% <https://spindynamics.org/wiki/index.php?title=diamond_frame_xyz.m>

function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)

% Check consistency
grumble(xaxis,yaxis,zaxis);

% Orthogonalise the supplied frame
frame=[xaxis(:) yaxis(:) zaxis(:)];
frame=diamond_frame_orth(frame);

end

% Consistency enforcement
function grumble(xaxis,yaxis,zaxis)
if(~isnumeric(xaxis)||~isreal(xaxis)||numel(xaxis)~=3)
    error('xaxis must be a real three-element array.');
end
if(~isnumeric(yaxis)||~isreal(yaxis)||numel(yaxis)~=3)
    error('yaxis must be a real three-element array.');
end
if(~isnumeric(zaxis)||~isreal(zaxis)||numel(zaxis)~=3)
    error('zaxis must be a real three-element array.');
end
end

% Experimental frames arrive approximate and leave orthogonal.

