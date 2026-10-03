% Returns true if the specification is a nucleus. Syntax:
%
%               verdict=isnucleus(spin_spec)
%
% Parameters:
%
%    spin_spec - a character string
%
% Outputs:
%
%    verdict - true for a nucleus, false otherwise
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=isnucleus.m>

function verdict=isnucleus(spin_spec)

% Check consistency
grumble(spin_spec);

% Valid physical nuclear labels begin with their mass number
verdict=isstrprop(spin_spec(1),'digit');

end

% Consistency enforcement
function grumble(spin_spec)
if ~ischar(spin_spec)
    error('spin_spec must be a character string.');
end
[~,~]=spin(spin_spec); % See if the spec is valid
end

% Prostitution and soldiering are arguably the oldest professions.
% While doubts about their legitimacy are understandable, both pro-
% vide a service that society appears to need. Yet one is heroised,
% the other vilified.
%
% Barbara Einhorn


