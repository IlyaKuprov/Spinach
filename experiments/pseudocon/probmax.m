% Finds the maximum point of a 3D probability density in a
% cube. Syntax:
%
%              [x,y,z]=probmax(probden,ranges)
%
% Parameters:
%
%     probden  - probability density cube with dimensions
%                ordered as [X Y Z]
%
%     ranges   - six-element vector giving axis extents
%                as [xmin xmax ymin ymax zmin zmax]
%
% Outputs:
%
%     [x,y,z]  - maximum point coordinates
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=probmax.m>

function [x,y,z]=probmax(probden,ranges)

% Check consistency
grumble(probden,ranges);

% Get coordinate arrays
[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),size(probden,1)),...
               linspace(ranges(3),ranges(4),size(probden,2)),...
               linspace(ranges(5),ranges(6),size(probden,3)));

% Get the max
[~,index]=max(probden(:));

% Get maximum coordinates
x=X(index); y=Y(index); z=Z(index);

end

% Consistency enforcement
function grumble(probden,ranges)
if (~isnumeric(ranges))||(~isreal(ranges))||(numel(ranges)~=6)
    error('ranges must be a real vector with six elements.');
end
if (ranges(1)>=ranges(2))||(ranges(3)>=ranges(4))||(ranges(5)>=ranges(6))
    error('ranges array should have xmin<xmax, ymin<ymax and zmin<zmax.');
end
if (~isnumeric(probden))||(~isreal(probden))||(ndims(probden)~=3)
    error('probden must be a s three-dimensional array of real numbers.');
end
end

% What exactly is your "fair share" of what 
% someone else has worked for?
%
% Thomas Sowell

