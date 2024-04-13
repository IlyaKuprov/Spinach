% Converts Cartesian coordinates [x y z] into spherical 
% coordinates according to the ISO convention. Syntax:
%
%              [r,theta,phi] = xyz2sph(x,y,z)
%
% Parameters:
%
%     x,y,z  - arrays of X, Y and Z coordinates
%
% Outputs:
%
%     r      - array of radii
%
%     theta  - array of inclinations
%
%     phi    - array of azimuth values
%
% e.suturina@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=xyz2sph.m>

function [r,theta,phi] = xyz2sph(x,y,z)

% Check consistency
grumble(x,y,z);

% Radius 0 <= r < Inf
r=sqrt(x.^2+y.^2+z.^2);

% Inclination 0 <= theta <= pi
theta=acos(z./r); 

% Azimuth 0 <= phi < 2*pi
phi=atan2(y,x);    

end

% Consistency enforcement
function grumble(x,y,z)
if (~isnumeric(x))||(~isnumeric(y))||(~isnumeric(z))
    error('all inputs must be numeric.');
end
if (~isreal(x))||(~isreal(y))||(~isreal(z))
    error('all inputs must be real.');
end
if (~all(size(x)==size(y)))||(~all(size(y)==size(z)))
    error('all inputs must have the same size.');
end
end

% Q: How would a theoretical physicist milk a cow?
% A: Consider a spherical cow in vacuum...
%
% a well-informed Christmas cracker

 