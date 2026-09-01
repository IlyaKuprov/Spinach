% Converts meV energy units used in solid state physics and
% phonon spectroscopy into Hz units preferred in magnetic
% resonance. Syntax:
%
%                       hz=mev2hz(mev)
%
% Arrays of any dimensions are supported. Parameters:
%
%   mev   - an array of values in milli-electronvolts
%
% Outputs:
%
%   hz    - an array of values in Hz
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=mev2hz.m>

function hz=mev2hz(mev)

% Check consistency
grumble(mev);

% Run the conversion
hz=1e-3*1.602176634e-19*mev/6.62607015e-34;

end

% Consistency enforcement
function grumble(mev)
if (~isnumeric(mev))||(~isreal(mev))
    error('the argument must be an array of real numbers.');
end
end

% O God, I could be bounded in a nutshell, and count
% myself a king of infinite space, were it not that I
% have bad dreams.
%
% William Shakespeare, Hamlet

