% Finds the extremum of a cubic interpolant built from function
% values and directional derivatives at two points and returns
% the best point inside the requested interpolation interval.
%
% Syntax:
%
%    [alpha,fx]=cubic_interp(end_A,end_B,alpha_A,alpha_B,...
%                            f_A,dir_deriv_A,f_B,dir_deriv_B)
%
% Arguments:
%
%    end_A             - first interpolation boundary in alpha space
%
%    end_B             - second interpolation boundary in alpha space
%
%    alpha_A           - first interpolation anchor point
%
%    alpha_B           - second interpolation anchor point
%
%    f_A               - function value at alpha_A
%
%    dir_deriv_A       - directional derivative at alpha_A
%
%    f_B               - function value at alpha_B
%
%    dir_deriv_B       - directional derivative at alpha_B
%
% Returns:
%
%    alpha             - selected maximiser of the cubic model
%
%    fx                - cubic model value at alpha
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cubic_interp.m>

function [alpha,fx]=cubic_interp(end_A,end_B,alpha_A,alpha_B,f_A,dir_deriv_A,f_B,dir_deriv_B)

% Check consistency
grumble(end_A,end_B,alpha_A,alpha_B,f_A,dir_deriv_A,f_B,dir_deriv_B);

% Build the cubic coefficients in normalised coordinates
c1=-2*(f_B-f_A)+(dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c2=3*(f_B-f_A)-(2*dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c3=(alpha_B-alpha_A)*dir_deriv_A; c4=f_A;

% Transform interpolation bounds to normalised coordinates
bounds=([end_A end_B]-alpha_A)/(alpha_B-alpha_A);

% Compute derivative roots of the cubic polynomial
s_points=roots([3*c1 2*c2 c3]);

% Remove complex roots from consideration
s_points(imag(s_points)~=0)=[];

% Remove roots lying outside interpolation bounds
s_points(s_points<min(bounds))=[];
s_points(s_points>max(bounds))=[];

% Include interpolation boundaries in candidate set
s_points=[min(bounds) s_points(:)' max(bounds)];

% Select the candidate with maximal cubic value
[fx,k]=max(polyval([c1 c2 c3 c4],s_points));

% Transform selected point back to alpha coordinates
alpha=alpha_A+s_points(k)*(alpha_B-alpha_A);

end

% Consistency enforcement
function grumble(end_A,end_B,alpha_A,alpha_B,f_A,dir_deriv_A,f_B,dir_deriv_B)
if isempty(end_A)||(~isnumeric(end_A))||(~isreal(end_A))||(~isscalar(end_A))
    error('end_A must be a real scalar.');
end
if isempty(end_B)||(~isnumeric(end_B))||(~isreal(end_B))||(~isscalar(end_B))
    error('end_B must be a real scalar.');
end
if isempty(alpha_A)||(~isnumeric(alpha_A))||(~isreal(alpha_A))||(~isscalar(alpha_A))
    error('alpha_A must be a real scalar.');
end
if isempty(alpha_B)||(~isnumeric(alpha_B))||(~isreal(alpha_B))||(~isscalar(alpha_B))
    error('alpha_B must be a real scalar.');
end
if alpha_A==alpha_B
    error('alpha_A and alpha_B must be different.');
end
if isempty(f_A)||(~isnumeric(f_A))||(~isreal(f_A))||(~isscalar(f_A))
    error('f_A must be a real scalar.');
end
if isempty(dir_deriv_A)||(~isnumeric(dir_deriv_A))||(~isreal(dir_deriv_A))||(~isscalar(dir_deriv_A))
    error('dir_deriv_A must be a real scalar.');
end
if isempty(f_B)||(~isnumeric(f_B))||(~isreal(f_B))||(~isscalar(f_B))
    error('f_B must be a real scalar.');
end
if isempty(dir_deriv_B)||(~isnumeric(dir_deriv_B))||(~isreal(dir_deriv_B))||(~isscalar(dir_deriv_B))
    error('dir_deriv_B must be a real scalar.');
end
end


