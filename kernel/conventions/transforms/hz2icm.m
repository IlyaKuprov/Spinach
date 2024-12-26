% Converts Hz units used in magnetic resonance into cm^-1 units 
% used in spectroscopy. Syntax:
%
%                         icm=hz2icm(hz)
%
% Arrays of any dimensions are supported. Parameters:
%
%   hz    - an array of values in Hz
%
% Outputs:
%
%   icm   - an array of values in inverse centimetres
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hz2icm.m>

function icm=hz2icm(hz)

% Check consistency
grumble(hz);

% Run the conversion
icm=hz/(100*299792458);

end

% Consistency enforcement
function grumble(hz)
if (~isnumeric(hz))||(~isreal(hz))
    error('the argument must be an array of real numbers.');
end
end

% "Damn, would I have to be nice to everyone for two years?!"
%
% IK, upon being informed that he was to
% organise the 48th ESR Group Conference
% in Southampton in two years' time.

