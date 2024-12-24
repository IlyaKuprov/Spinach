% Generates a convex hull of a two-angle grid for 2D 
% surface plotting. Syntax:
%
%       hull=get_hull(theta_angles,phi_angles)
%
% Parameters:
%
%   theta_angles - a column vector of theta angles,
%                  polar coordinates, ISO convention,
%                  radians
%
%   phi_angles   - a column vector of phi angles,
%                  polar coordinates, ISO convention,
%                  radians
%
% Outputs:
%
%   hull         - a matrix of point indices of 
%                  dimension Nx3, where N is the
%                  number of triangular facets
%
%   edges        - a matrix of point indices of 
%                  dimension Nx2, where N is the
%                  number of grid edges
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=get_hull.m>

function [hull,edges]=get_hull(theta_angles,phi_angles)

% Check consistency
grumble(theta_angles,phi_angles)

% Get Cartesian coordinates
x=sin(theta_angles).*cos(phi_angles);
y=sin(theta_angles).*sin(phi_angles);
z=cos(theta_angles);

% Get the convex hull
hull=convhull(x,y,z);

% Get the edges
edges=unique([hull(:,1) hull(:,2);
              hull(:,2) hull(:,3)],'rows');
edges=unique([edges; edges(:,[2 1])],'rows');
edges(edges(:,1)==edges(:,2),:)=[];

end

% Consistency enforcement
function grumble(theta_angles,phi_angles)
if (~isnumeric(theta_angles))||(~isreal(theta_angles))||...
   any(~isfinite(theta_angles))||(size(theta_angles,2)~=1)
    error('theta_angles must be a column vector of real numbers.');
end
if (~isnumeric(phi_angles))||(~isreal(phi_angles))||...
   any(~isfinite(phi_angles))||(size(phi_angles,2)~=1)
    error('phi_angles must be a column vector of real numbers.');
end
if numel(theta_angles)~=numel(phi_angles)
    error('the number of elements in theta_angles and phi_angles must be the same.');
end
end

% But let me offer you my definition of social justice: I keep
% what I earn and you keep what you earn. Do you disagree? Well
% then tell me how much of what I earn belongs to you - and why?
%
% Walter E. Williams, "All It Takes is Guts: a Minority View"

