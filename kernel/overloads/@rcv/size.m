% Returns the size of an RCV object. Syntax:
%
%                       s=size(obj)
%                       s=size(obj,dim)
%
% Parameters:
%
%    obj   - RCV object
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

% Return either both dimensions or a specific one
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
    error('the first argument must be an rcv object.');
end
if ~isscalar(obj)
    error('the rcv input must be a scalar object.');
end
if nargin<2
    dim=[];
end
if (nargin==2)&&~(isscalar(dim)&&isnumeric(dim))
    error('dimension index must be a numeric scalar.');
end
end

% I did not succeed in life by intelligence. I succeeded
% because I have a long attention span.
%
% Charlie Munger
