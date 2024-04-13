% Arc length between two points on the unit sphere specified by
% the unit vectors supplied. Syntax:
%
%                     sig=arclength(r1,r2)
%
% Parameters:
%
%   r1,r2    - three-element unit vectors with Cartesian
%              coordinates of the arc endpoints
%
% Outputs:
%
%   sig      - arc length
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=arclength.m>

function sig=arclength(r1,r2)

% Check consistency
grumble(r1,r2);

% Normalise the vectors
r1=r1/norm(r1,2); r1=r1(:);
r2=r2/norm(r2,2); r2=r2(:);

% Get the arc length
sig=atan2(norm(cross(r1,r2),2),dot(r1,r2));

end

% Consistency enforcement
function grumble(r1,r2)
if (~isnumeric(r1))||(~isreal(r1))||(numel(r1)~=3)||...
   (~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)
    error('r1,r2 must be three-element real vectors.');
end
if (abs(norm(r1,2)-1)>sqrt(eps))||...
   (abs(norm(r2,2)-1)>sqrt(eps))
    error('r1,r2 must be unit vectors.');
end
end

% Ah, there's nothing more exciting than science. You get all the
% fun of sitting still, being quiet, writing down numbers, paying
% attention... science has it all.
%
% Principal Skinner

