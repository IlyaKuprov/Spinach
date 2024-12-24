% Converts cm^-1 units used in spectroscopy into Hz units 
% preferred in magnetic resonance. Syntax:
%
%                     hz=icm2hz(icm)
%
% Arrays of any dimensions are supported. Parameters:
%
%   icm   - an array of values in inverse centimetres
%
% Outputs:
%
%   hz    - an array of values in Hz
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=icm2hz.m>

function hz=icm2hz(icm)

% Check consistency
grumble(icm);

% Run the conversion
hz=100*299792458*icm;

end

% Consistency enforcement
function grumble(icm)
if (~isnumeric(icm))||(~isreal(icm))
    error('the argument must be an array of real numbers.');
end
end

% Satire thrives where the usual checks on human 
% folly fail.
%
% The Economist

