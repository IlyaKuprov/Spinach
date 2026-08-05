% Converts Kelvin energy units used for Debye temperatures and
% thermal energy scales in solid state physics into Hz units
% preferred in magnetic resonance. Syntax:
%
%                     hz=kelvin2hz(kelvin)
%
% Arrays of any dimensions are supported. Parameters:
%
%   kelvin  - an array of values in Kelvin
%
% Outputs:
%
%   hz      - an array of values in Hz
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kelvin2hz.m>

function hz=kelvin2hz(kelvin)

% Check consistency
grumble(kelvin);

% Run the conversion
hz=1.380649e-23*kelvin/6.62607015e-34;

end

% Consistency enforcement
function grumble(kelvin)
if (~isnumeric(kelvin))||(~isreal(kelvin))
    error('the argument must be an array of real numbers.');
end
end

% When you can measure what you are speaking about, and
% express it in numbers, you know something about it; but
% when you cannot measure it, when you cannot express it
% in numbers, your knowledge is of a meagre and unsatis-
% factory kind.
%
% William Thomson, Lord Kelvin


