% Returns true if the particle specification is not an electron,
% a neutron, a muon, or a ghost. Syntax:
%
%                  verdict=isnucleus(spin_spec)
%
% Parameters:
%
%    spin_spec - a character string
%
% Output:
%
%    verdict - true for a nucleus known to Spinach, false
%              for an electron, a muon, a neutron, or a ghost
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=isnucleus.m>

function verdict=isnucleus(spin_spec)

% Check consistency
grumble(spin_spec);

% A simple matching check
if strcmp(spin_spec(1),'E')||...
   ismember(spin_spec,{'G','E','N','M'})
    verdict=false();
else
    verdict=true();
end

end

% Consistency enforcement
function grumble(spin_spec)
if ~ischar(spin_spec)
    error('spin_spec must be a character string.');
end
[~,~]=spin(spin_spec);
end

% Prostitution and soldiering are arguably the oldest professions.
% While doubts about their legitimacy are understandable, both pro-
% vide a service that society appears to need. Yet one is heroised,
% the other vilified.
%
% Barbara Einhorn

