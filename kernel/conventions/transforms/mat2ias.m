% Isotropic-antisymmetric-symmetric decomposition of a 3x3
% real interaction matrix between real vectors u and v:
%
%        u'*C*v = a*(u'*v) + d'*cross(u,v) + u'*A*v
% 
% Syntax:
%
%                      [a,d,A]=mat2ias(C)
%
% Parameters:
%
%    C - real 3x3 matrix
%
% Outputs:
%
%    a - scalar component
%
%    d - antisymmetric coupling vector
%
%    A - symmetric coupling matrix
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mat2ias.m>

function [a,d,A]=mat2ias(C)

% Check consistency
grumble(C);

% Isotropic part
a=trace(C)/3;

% Antisymmetric part
d=[C(2,3)-C(3,2);
   C(3,1)-C(1,3);
   C(1,2)-C(2,1)]/2;

% Traceless symmetric part
A=(C+C')/2-a*eye(3,3);

end

% Consistency enforcement
function grumble(C)
if (~isnumeric(C))||(~isreal(C))||...
   (size(C,1)~=3)||(size(C,2)~=3)
    error('C must be a real 3x3 matrix.');
end
end

% Mathematics is the only true metaphysics. 
%
% Lord Kelvin

