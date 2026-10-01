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
grumble(x,dim);

% Implicit core descriptions carry their dimensions
if isstruct(x)
    n=x.dims(dim);
else
    n=size(x,dim);
end

end

% Consistency enforcement
function grumble(x,dim)
if isstruct(x)
    if ~isscalar(x)||~isfield(x,'dims')
        error('implicit core descriptions need a dims field.');
    end
    if ~isnumeric(x.dims)||~isreal(x.dims)||~isequal(size(x.dims),[1 2])||...
       any(~isfinite(x.dims))||any(x.dims<1)||any(mod(x.dims,1)~=0)
        error('implicit core dims must be a row of two positive integers.');
    end
elseif ~isnumeric(x)||~ismatrix(x)
    error('x must be a matrix or an implicit core description.');
end
if (~isnumeric(dim))||(~isscalar(dim))||(~ismember(dim,[1 2]))
    error('dim must be 1 or 2.');
end
end

% Simplicity is the ultimate sophistication.
%
% Leonardo da Vinci


