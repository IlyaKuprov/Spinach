% Solid angle of a convex spherical polygon as described in
%  
%         https://doi.org/10.1109/TBME.1983.325207
%
% Syntax:
%
%              A=one_vcell_solidangle(v,centre)
%
% Parameters:
%
%     v      - (3 x n) matrix of unit vectors giving
%               the coordinates of each vertex
%
%     centre - centre vertex coordinates, optional
%
% Outputs:
%
%     S      - the solid angle, radians
%
% brunoluong@yahoo.com
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=one_vcell_solidangle.m>

function S=one_vcell_solidangle(v,centre)

% Check consistency
grumble(v,centre);

% Straightforward math
if nargin<2
    n=size(v,2);
    s=zeros(1,n-2);
    for k=2:n-1
        T=v(:,[1 k k+1]);
        num=det(T);
        denom=1+sum(sum(T.*T(:,[2 3 1]),1),2);
        s(k-1)=num/denom;
    end
else
    v(:,end+1)=v(:,1);
    n=size(v,2);
    s=zeros(1,n-1);
    for k=1:n-1
        T=[centre, v(:,[k k+1])];
        num=det(T);
        denom=1+sum(sum(T.*T(:,[2 3 1]),1),2);
        s(k)=num/denom;
    end
end
S=atan(s); S=2*sum(S);

end

% Consistency enforcement
function grumble(v,centre)
if any(abs(sum(v.^2,1)-1)>1e-6)
    error('v must contain unit vectors.');
end
if abs(norm(centre,2)-1)>1e-6
    error('centre muslt be a unit vector.');
end
end

% Being a mathematician is a bit like being a manic 
% depressive: you spend your life alternating between
% giddy elation and black despair.
%
% Steven G. Krantz

