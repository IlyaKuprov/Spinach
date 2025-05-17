% Returns true if the particle is an electron. Syntax:
%
%             verdict=iselectron(spin_spec)
%
% Parameters:
%
%    spin_spec - a Spinach particle specification
%
% Output:
%
%    verdict - true for an electron, false otherwise
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=iselectron.m>

function verdict=iselectron(spin_spec)

% Check consistency
grumble(spin_spec);

% A simple matching check
if ismember(spin_spec(1),{'E'})
    verdict=true();
else
    verdict=false();
end

end

% Consistency enforcement
function grumble(spin_spec)
if ~ischar(spin_spec)
    error('spin_spec must be a character string.');
end
[~,~]=spin(spin_spec); % See if the spec is valid
end

% Лучше умереть героем чем жить пидорасом.
%
% Евгений Пригожин

