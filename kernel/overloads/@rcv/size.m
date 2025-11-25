% Returns the size of an RCV sparse matrix. Syntax:
%
%                    s=size(A,dim)
%                      
%
% Parameters:
%
%    A     - RCV sparse matrix
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

function s=size(A,dim)

% Check consistency
if nargin==1
    grumble(A);
else
    grumble(A,dim);
end

% Mimic Matlab
if nargin==1
    s=[A.numRows A.numCols];
else
    if dim==1
        s=A.numRows;
    elseif dim==2
        s=A.numCols;
    else
        s=1;
    end
end

end

% Consistency enforcement
function grumble(A,dim)
if ~isa(A,'rcv')
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

