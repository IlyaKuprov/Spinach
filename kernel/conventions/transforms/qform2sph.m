% Returns the spherical harmonic expansion coefficients of the
% following quadratic form:
%
%    [x y z]*A*[x y z]'/norm([x y z],2)^2 = sum(r_LM*Y_LM)
%
% Syntax:
%
%                   [r0,r1,r2]=qform2sph(A)
% 
% Parameters:
%
%   A - a symmetric 3x3 matrix
%
% Output:
% 
%   [r0,r1,r2] - coefficients for zero, first, 
%                and second rank spherical
%                harmonics in the order of
%                decreasing m index of Ylm
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=qform2sph.m>

function [r0,r1,r2]=qform2sph(A)

% Check consistency
grumble(A);

% Zeroth rank
r0=(2/3)*sqrt(pi)*trace(A);

% First rank from 1 to -1
r1(1)=0; r1(2)=0; r1(3)=0;

% Second rank from 2 to -2
r2(1)=+sqrt(2*pi/15)*(A(1,1)-A(2,2)-1i*(A(1,2)+A(2,1)));
r2(2)=-sqrt(2*pi/15)*(A(1,3)+A(3,1)-1i*(A(2,3)+A(3,2)));
r2(3)=-(2/3)*sqrt(pi/5)*(-2*A(3,3)+A(2,2)+A(1,1));
r2(4)=+sqrt(2*pi/15)*(A(1,3)+A(3,1)+1i*(A(2,3)+A(3,2)));
r2(5)=+sqrt(2*pi/15)*(A(1,1)-A(2,2)+1i*(A(1,2)+A(2,1)));

end

% Consistency enforcement
function grumble(A)
if (~isnumeric(A))||(size(A,1)~=3)||(size(A,2)~=3)||...
   (~isreal(A))||(~issymmetric(A))
    error('A must be a real symmetric 3x3 matrix.');
end
end

% Boring things in physics may be defined as 
% those that make sense.
%
% Source unknown

