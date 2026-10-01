% Dimension of a polyadic core, which may be a matrix, a polyadic,
% an opium, or an implicit core description with a dims field. Syntax:
%
%                         n=core_size(x,dim)
%
% Parameters:
%
%    x   - a polyadic core
%
%    dim - dimension whose size is required, 1 or 2
%
% Outputs:
%
%    n   - the number of rows or columns of the core
%
% ilya.kuprov@weizmann.ac.il

function n=core_size(x,dim)

% Check consistency
grumble(dim);

% Implicit core descriptions carry their dimensions
if isstruct(x)
    n=x.dims(dim);
else
    n=size(x,dim);
end

end

% Consistency enforcement
function grumble(dim)
if (~isscalar(dim))||(~ismember(dim,[1 2]))
    error('dim must be 1 or 2.');
end
end

% Simplicity is the ultimate sophistication.
%
% Leonardo da Vinci


