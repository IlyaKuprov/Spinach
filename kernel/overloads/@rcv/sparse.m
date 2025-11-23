% Converts an RCV sparse matrix into a Matlab sparse 
% matrix. Syntax:
%
%                     S=sparse(obj)
%
% Parameters:
%
%    obj   - RCV sparse matrix
%
% Outputs:
%
%    S     - Matlab sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/sparse.m>

function S=sparse(obj)

% Check consistency
grumble(obj);

% Check if empty
if isempty(obj.col)

    % Empty matrix of a specified size
    S=spalloc(obj.numRows,obj.numCols,0);

else

    % Call Matlab's sparse matrix constructor
    S=sparse(obj.row,obj.col,obj.val,obj.numRows,obj.numCols);

end

end

% Consistency enforcement
function grumble(obj)
if ~isa(obj,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% Working 16 hours a day, 7 days a week, 52 weeks
% in a year, and people still calling me lucky.
%
% Elon Musk

