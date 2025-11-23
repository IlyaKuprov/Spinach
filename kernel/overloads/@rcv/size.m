% Returns the size of an RCV sparse matrix. Syntax:
%
%                    s=size(obj,dim)
%                      
%
% Parameters:
%
%    obj   - RCV object
%
%    dim   - optional dimension index
%
% Outputs:
%
%    s     - size vector or dimension length
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/size.m>

function s=size(obj,dim)

% Check consistency
if nargin==1
    grumble(obj);
else
    grumble(obj,dim);
end

% Mimic Matlab
if nargin==1
    s=[obj.numRows obj.numCols];
else
    if dim==1
        s=obj.numRows;
    elseif dim==2
        s=obj.numCols;
    else
        s=1;
    end
end

end

% Consistency enforcement
function grumble(obj,dim)
if ~isa(obj,'rcv')
    error('the first argument must be an RCV sparse matrix.');
end
if nargin==2
    if (~isscalar(dim))||(~isnumeric(dim))
        error('dimension index must be a numeric scalar.');
    end
end
end

% I did not succeed in life by intelligence. I succeeded
% because I have a long attention span.
%
% Charlie Munger

