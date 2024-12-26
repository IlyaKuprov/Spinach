% Converts resonance offsets into chemical shifts. Syntax:
%
%                   ppm=hz2ppm(hz,B0,nucleus)
%
% Inputs:
%
%    hz        - resonance offset in Hz
%
%    B0        - magnet induction, Tesla
%
%    nucleus   - a string specifying the isotope, e.g. '1H'
%
% Result:
%
%    ppm       - chemical shift in ppm
%
% Note: signs of the magnetogyric ratios are preserved.
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hz2ppm.m>

function ppm=hz2ppm(hz,B0,nucleus)

% Check consistency
grumble(hz,B0,nucleus);

% Calculate chemical shift in Hz
ppm=1e6*(2*pi*hz)/(B0*spin(nucleus)); 

end

% Consistency enforcement
function grumble(hz,B0,nucleus)
if (~isnumeric(hz))||(~isreal(hz))
    error('resonance offset must be real.');
end
if (~isnumeric(B0))||(~isreal(B0))
    error('magnetic induction must be real.');
end
if ~ischar(nucleus)
    error('nucleus must be a character array');
end
end

% Somebody Else's Wife: oh, I am so tired of being boringly taken
%                       for granted, he never really acknowledges
%                       me anymore, I want to have a turbulent af-
%                       fair with a poet!
%
% IK: sorry, we are out of poets - would a scientist do?
%
% Somebody Else's Wife: yes... 

