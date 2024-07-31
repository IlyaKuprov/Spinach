% Multiplies all entries of a structure by a user-specified mat-
% rix. Nested structures are processed recursively. Syntax:
%
%                     str_out=mtimes(M,str_in)
%
% Parameters:
%
%     M - any numeric object (scalar, matrix, etc.)
%
%     str_in - a structure with numeric subfields
%
% Outputs:
%
%     str_out - the resulting structure
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=struct/mtimes.m>

function str_out=mtimes(M,str_in)

% Check consistency
grumble(M,str_in)

% Get the field names
fnames=fieldnames(str_in);

% Loop over field names
for n=1:numel(fnames)
    
    % Recursive call for each field name
    str_out.(fnames{n})=M*str_in.(fnames{n});
    
end

end

% Consistency enforcement
function grumble(M,str_in)
if ~isnumeric(M)
    error('the first argument must be numeric.');
end
if ~isstruct(str_in)
    error('the second argument must be a structure.');
end
end

% Arthur Dent:  What happens if I press this button?
% Ford Prefect: I wouldn't --
% Arthur Dent:  Oh.
% Ford Prefect: What happened?
% Arthur Dent:  A sign lit up, saying "please do not press this button again".
%
% Douglas Adams

