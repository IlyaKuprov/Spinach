% Draws a horizontal stripe of 10% grey (black at 10%
% opacity) between the specified Y coordinates in the
% current axis system; a shaded strip analogue of the
% yline function in Matlab. The stripe spans the full
% width of the axes and sits underneath the data.
% Syntax:
%
%                     khorzshade(a,b)
%
% Parameters:
%
%    a, b - real numbers, Y coordinates of the
%           two edges of the stripe, in the units
%           of the Y axis
%
% Outputs:
%
%    updates the current axis system
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=khorzshade.m>

function khorzshade(a,b)

% Check consistency
grumble(a,b);

% Black stripe at 10% opacity
yregion(a,b,'FaceColor',[0 0 0],'FaceAlpha',0.1);

end

% Consistency enforcement
function grumble(a,b)
if (~isnumeric(a))||(~isreal(a))||(~isscalar(a))||(~isfinite(a))
    error('a must be a finite real number.');
end
if (~isnumeric(b))||(~isreal(b))||(~isscalar(b))||(~isfinite(b))
    error('b must be a finite real number.');
end
end

