% Converts an RCV object to a sparse matrix. Syntax:
%
%                       S=sparse(obj)
%
% Parameters:
%
%    obj   - RCV object
%
% Outputs:
%
%    S     - sparse matrix representation of obj
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/sparse.m>

function S=sparse(obj)

% Check consistency
grumble(obj);

% Build a sparse matrix from row, column and value arrays
if isempty(obj.col)
    S=spalloc(obj.numRows,obj.numCols,0);
else
    S=sparse(obj.row,obj.col,obj.val,obj.numRows,obj.numCols);
end

end

% Consistency enforcement
function grumble(obj)
if ~isa(obj,'rcv')
    error('the input must be an rcv object.');
end
if ~isscalar(obj)
    error('the input must be a scalar rcv object.');
end
end

% Working 16 hours a day, 7 days a week, 52 weeks
% in a year, and people still calling me lucky.
%
% Elon Musk
