% Field getter with a default value. Syntax:
%
%          value=diamond_get(parameters,name,default)
%
% This internal helper returns a structure field when present and a
% specified default otherwise.
%
% Parameters:
%
%    parameters  - input structure
%
%    name        - field name
%
%    default     - default value
%
% Outputs:
%
%    value       - field value or default
%
% <https://spindynamics.org/wiki/index.php?title=diamond_get.m>

function value=diamond_get(parameters,name,default)

% Check consistency
grumble(parameters,name);

% Return the field or its default
if isfield(parameters,name)
    value=parameters.(name);
else
    value=default;
end

end

% Consistency enforcement
function grumble(parameters,name)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if ~ischar(name)
    error('name must be a character string.');
end
end

% Defaults are useful only when they are explicit.

