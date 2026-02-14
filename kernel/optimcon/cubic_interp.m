% Finds the extremum of a cubic interpolant built from function
% values and directional derivatives at two points and returns
% the best point inside the interpolation interval. Syntax:
%
%   [alpha,fx]=cubic_interp(end_A,end_B,alpha_A,alpha_B,...
%                           f_A,dir_deriv_A,f_B,dir_deriv_B)
%
% Arguments:
%
%    end_a      - first interpolation boundary in alpha space
%
%    end_b      - second interpolation boundary in alpha space
%
%    alpha_a    - first interpolation anchor point
%
%    alpha_b    - second interpolation anchor point
%
%    f_a        - function value at alpha_a
%
%    dir_der_a  - directional derivative at alpha_a
%
%    f_b        - function value at alpha_b
%
%    dir_der_b  - directional derivative at alpha_b
%
% Returns:
%
%    alpha      - selected maximiser of the cubic model
%
%    fx         - cubic model value at alpha
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cubic_interp.m>

function [alpha,fx]=cubic_interp(end_a,end_b,alpha_a,alpha_b,...
                                 f_a,dir_der_a,f_b,dir_der_b)

% Check consistency
grumble(end_a,end_b,alpha_a,alpha_b,f_a,dir_der_a,f_b,dir_der_b);

% Build the cubic coefficients in normalised coordinates
c1=-2*(f_b-f_a)+(dir_der_a+dir_der_b)*(alpha_b-alpha_a);
c2=3*(f_b-f_a)-(2*dir_der_a+dir_der_b)*(alpha_b-alpha_a);
c3=(alpha_b-alpha_a)*dir_der_a; c4=f_a;

% Transform interpolation bounds to normalised coordinates
bounds=([end_a end_b]-alpha_a)/(alpha_b-alpha_a);

% Compute derivative roots 
s_points=roots([3*c1 2*c2 c3]);

% Remove complex roots 
s_points(imag(s_points)~=0)=[];

% Remove roots outside interpolation bounds
s_points(s_points<min(bounds))=[];
s_points(s_points>max(bounds))=[];

% Include interpolation boundaries in candidate set
s_points=[min(bounds) s_points(:)' max(bounds)];

% Select the candidate with maximal cubic value
[fx,k]=max(polyval([c1 c2 c3 c4],s_points));

% Transform selected point back to alpha coordinates
alpha=alpha_a+s_points(k)*(alpha_b-alpha_a);

end

% Consistency enforcement
function grumble(end_a,end_b,alpha_a,alpha_b,f_a,dir_der_a,f_b,dir_der_b)
if isempty(end_a)||(~isnumeric(end_a))||(~isreal(end_a))||(~isscalar(end_a))
    error('end_A must be a real scalar.');
end
if isempty(end_b)||(~isnumeric(end_b))||(~isreal(end_b))||(~isscalar(end_b))
    error('end_B must be a real scalar.');
end
if isempty(alpha_a)||(~isnumeric(alpha_a))||(~isreal(alpha_a))||(~isscalar(alpha_a))
    error('alpha_A must be a real scalar.');
end
if isempty(alpha_b)||(~isnumeric(alpha_b))||(~isreal(alpha_b))||(~isscalar(alpha_b))
    error('alpha_B must be a real scalar.');
end
if alpha_a==alpha_b
    error('alpha_A and alpha_B must be different.');
end
if isempty(f_a)||(~isnumeric(f_a))||(~isreal(f_a))||(~isscalar(f_a))
    error('f_A must be a real scalar.');
end
if isempty(dir_der_a)||(~isnumeric(dir_der_a))||(~isreal(dir_der_a))||(~isscalar(dir_der_a))
    error('dir_deriv_A must be a real scalar.');
end
if isempty(f_b)||(~isnumeric(f_b))||(~isreal(f_b))||(~isscalar(f_b))
    error('f_B must be a real scalar.');
end
if isempty(dir_der_b)||(~isnumeric(dir_der_b))||(~isreal(dir_der_b))||(~isscalar(dir_der_b))
    error('dir_deriv_B must be a real scalar.');
end
end

% It is baked beans that should be banned.
% In nauseating sauce encanned,
% All sickly-sweet and liquidous
% The bastards are ubiquitous.
%
% An English Breakfast should eschew
% This transatlantic parvenu –
% No beans! No beans! I cry. Too late.
% Plop goes the dollop on my plate.
%
% Not even a mixed grill avoids
% The devil's horrid haemorrhoids.
% I scrape them to the side, but still
% They manage to pollute my meal.
%
% They sidle up and touch my chips,
% Contriving thus to reach my lips
% So I again am forced to taste
% Their overtones of toxic waste.
%
% Ann Drysdale

