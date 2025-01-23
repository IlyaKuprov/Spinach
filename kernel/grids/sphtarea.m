% Area of the curvilinear triangle on the unit sphere defined 
% by the vertex coordinates supplied. Syntax:
%
%                  S=sphtarea(r1,r2,r3,sflag)
%
% Parameters:
%
%   r1,r2,r3  - three-element unit vectors with Cartesian
%               coordinates of triangle vertices
%
%      sflag  - 'signed' would take into account surface
%               normal direction, 'unsigned' (default)
%               would always return a positive area
%
% Outputs:
%
%          S  - spherical triangle surface area
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sphtarea.m>

function S=sphtarea(r1,r2,r3,sflag)

% Default to unsigned area
if ~exist('sflag','var'), sflag='unsigned'; end

% Check consistency
grumble(r1,r2,r3,sflag);

% Stretch the vectors
r1=r1(:); r2=r2(:); r3=r3(:);

% Get signed area
S=2*atan2(det([r1 r2 r3]),dot(r1,r2)+dot(r2,r3)+dot(r3,r1)+1);

% Kill the sign if appropriate
if strcmp(sflag,'unsigned'), S=abs(S); end

end

% Consistency enforcement
function grumble(r1,r2,r3,sflag)
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
if (arclength(r1,r2)>pi/2)||...
   (arclength(r2,r3)>pi/2)||...
   (arclength(r3,r1)>pi/2)
    error('triangle arc lengths cannot exceed pi/2.');
end
if (~ischar(sflag))||...
   (~ismember(sflag,{'signed','unsigned'}))
    error('sflag must be ''signed'' or ''unsigned''.');
end
end

% Коллега, глядя на испорченный лаптоп: "Сказать что кота вырвало,
% это значит снять ответственность с кота и возложить ее на какие-
% то неведомые силы. Кота не вырвало, кот - наблевал!"

