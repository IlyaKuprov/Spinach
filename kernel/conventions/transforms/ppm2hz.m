% Converts chemical shifts into resonance offsets. Syntax:
%
%                   hz=ppm2hz(ppm,B0,nucleus)
%
% Inputs:
%
%    ppm       - chemical shift in ppm
%
%    B0        - magnet induction, Tesla
%
%    nucleus   - a string specifying the isotope, e.g. '1H'
%
% Result:
%
%    hz        - resonance offset in Hz
%
% Note: signs of the magnetogyric ratios are preserved.
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ppm2hz.m>

function hz=ppm2hz(ppm,B0,nucleus)

% Check consistency
grumble(ppm,B0,nucleus);

% Calculate chemical shift in Hz
hz=1e-6*ppm*(B0*spin(nucleus)/(2*pi)); 

end

% Consistency enforcement
function grumble(ppm,B0,nucleus)
if (~isnumeric(ppm))||(~isreal(ppm))
    error('chemical shift must be real.');
end
if (~isnumeric(B0))||(~isreal(B0))
    error('magnetic induction must be real.');
end
if ~ischar(nucleus)
    error('nucleus must be a character array');
end
end

% "Smoking - NO HYDROGEN!"
%
% Safety warning on Anatole Abragam's door

