% Converts magnetic susceptibility from the Angstrom^3 units 
% required by Spinach pseudocontact shift functionality into
% the cgs-ppm (aka cm^3/mol) units quoted by quantum chemist-
% ry packages. Syntax:
%
%                  cgsppm=ang2cgsppm(ang)
%
% Parameters:
%
%   ang    - an array of values in cubic Angstrom
%
% Outputs:
%
%   cgsppm - an array of values in cgs-ppm
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ang2cgsppm.m>

function cgsppm=ang2cgsppm(ang)

% Check consistency
grumble(ang);

% Do the calculation
cgsppm=6.02214129e23*ang/(4*pi*1e18);

end

% Consistency enforcement
function grumble(ang)
if ~isnumeric(ang)
    error('input must be numeric.');
end
end

% No artist tolerates reality.
%
% Friedrich Nietzsche

