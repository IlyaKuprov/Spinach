% Unit vector from polar angles in the diamond frame. Syntax:
%
%          v=diamond_sph_vec(theta,phi)
%
% This internal helper returns a unit vector for polar angles specified
% in degrees.
%
% Parameters:
%
%    theta  - polar angle, degrees
%
%    phi    - azimuthal angle, degrees
%
% Outputs:
%
%    v      - unit vector
%
% <https://spindynamics.org/wiki/index.php?title=diamond_sph_vec.m>

function v=diamond_sph_vec(theta,phi)

% Check consistency
grumble(theta,phi);

% Convert spherical angles to a Cartesian unit vector
v=[sind(theta)*cosd(phi);sind(theta)*sind(phi);cosd(theta)];

end

% Consistency enforcement
function grumble(theta,phi)
if(~isnumeric(theta)||~isreal(theta)||~isscalar(theta))
    error('theta must be a real scalar.');
end
if(~isnumeric(phi)||~isreal(phi)||~isscalar(phi))
    error('phi must be a real scalar.');
end
end

% Angles are data only when their convention is explicit.

