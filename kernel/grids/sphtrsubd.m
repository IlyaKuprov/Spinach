% Spherical triangle subdivision. Returns the midpoints of 
% the sides of a spherical triangle specified by the unit
% vectors supplied. Syntax:
%
%             [r12,r23,r31]=sphtrsubd(r1,r2,r3)
%
% Parameters:
%
%   r1,r2,r3    - three-element unit vectors with Cartesian
%                 coordinates of triangle vertices
%
% Outputs:
%
%   r12,r23,r31 - three-element unit vectors with Cartesian
%                 coordinates of triangle arc midpoints
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sphtrsubd.m>

function [r12,r23,r31]=sphtrsubd(r1,r2,r3)

% Check consistency
grumble(r1,r2,r3);

% Not particularly hard
r12=r1+r2; r12=r12/norm(r12,2);
r23=r2+r3; r23=r23/norm(r23,2);
r31=r3+r1; r31=r31/norm(r31,2);

end

% Consistency enforcement
function grumble(r1,r2,r3)
if (~isnumeric(r1))||(~isreal(r1))||(numel(r1)~=3)||...
   (~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)||...
   (~isnumeric(r3))||(~isreal(r3))||(numel(r3)~=3)
    error('r1,r2,r3 must be three-element real vectors.');
end
if (abs(norm(r1,2)-1)>sqrt(eps))||...
   (abs(norm(r2,2)-1)>sqrt(eps))||...
   (abs(norm(r3,2)-1)>sqrt(eps))
    error('r1,r2,r3 must be unit vectors.');
end
if sphtarea(r1,r2,r3)>pi/2
    error('parent triangle area cannot exceed pi/2.');
end
if (arclength(r1,r2)>pi/2)||...
   (arclength(r2,r3)>pi/2)||...
   (arclength(r3,r1)>pi/2)
    error('parent triangle arc lengths cannot exceed pi/2.');
end
end

% Feminism was established so as to allow unattractive women
% easier access to the mainstream society.
%
% Rush Limbaugh

