% Zooms a 3D data cube to the fractional limits specified 
% by the user. Syntax:
%
%      [density,ext]=zoom_3d(density,ext,zoom_ranges)
%
% Parameters:
%
%      density  - probability density cube with dimensions
%                 ordered as [X Y Z]
%
%      ext      - grid extents in Angstrom, ordered as
%                 [xmin xmax ymin ymax zmin zmax]
%
%   zoom_ranges - zoom ranges along each axis as fractions,
%                 ordered as [xmin xmax ymin ymax zmin zmax],
%                 e.g. [0.3 0.6 0.1 0.2 0.5 0.8]
%
% Outputs:
%
%      density  - probability density cube with dimensions
%                 ordered as [X Y Z]
%
%      ext      - grid extents in Angstrom, ordered as
%                 [xmin xmax ymin ymax zmin zmax]
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=zoom_3d.m>

function [density,ext]=zoom_3d(density,ext,zoom_ranges)

% Check consistency
grumble(density,ext,zoom_ranges);

% Generate axis ticks
oldxvals=linspace(ext(1),ext(2),size(density,1));
oldyvals=linspace(ext(3),ext(4),size(density,2));
oldzvals=linspace(ext(5),ext(6),size(density,3));

% Get new range indices
xmin=max([1 floor(size(density,1)*zoom_ranges(1))]);
ymin=max([1 floor(size(density,2)*zoom_ranges(3))]);
zmin=max([1 floor(size(density,3)*zoom_ranges(5))]);
xmax=min([size(density,1) ceil(size(density,1)*zoom_ranges(2))]);
ymax=min([size(density,2) ceil(size(density,2)*zoom_ranges(4))]);
zmax=min([size(density,3) ceil(size(density,3)*zoom_ranges(6))]);

% Extract the subcube
density=density(xmin:xmax,ymin:ymax,zmin:zmax);

% Update the extents
ext=[oldxvals(xmin) oldxvals(xmax) ...
     oldyvals(ymin) oldyvals(ymax) ...
     oldzvals(zmin) oldzvals(zmax)];

end

% Consistency enforcement
function grumble(density,ext,zoom_ranges)
if (~isnumeric(ext))||(~isreal(ext))||(numel(ext)~=6)
    error('ext must be a real vector with six elements.');
end
if (ext(1)>=ext(2))||(ext(3)>=ext(4))||(ext(5)>=ext(6))
    error('ext array should have xmin<xmax, ymin<ymax and zmin<zmax.');
end
if (~isnumeric(density))||(~isreal(density))||(ndims(density)~=3)
    error('density must be a three-dimensional array of real numbers.');
end
if (~isnumeric(zoom_ranges))||(~isreal(zoom_ranges))||(numel(zoom_ranges)~=6)
    error('zoom_ranges must be a real vector with six elements.');
end
if (zoom_ranges(1)>=zoom_ranges(2))||...
   (zoom_ranges(3)>=zoom_ranges(4))||...
   (zoom_ranges(5)>=zoom_ranges(6))
    error('zoom_ranges array should have xmin<xmax, ymin<ymax and zmin<zmax.');
end
if any(zoom_ranges<0)||any(zoom_ranges>1)
    error('all elements of zoom_ranges must be between 0 and 1.');
end
end

% Either you think - or else others have to think for you and
% take power from you, pervert and discipline your natural ta-
% stes, civilize and sterilize you.
%
% F. Scott Fitzgerald

