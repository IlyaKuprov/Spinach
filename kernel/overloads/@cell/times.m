% Multiplies all entries of a cell array by a user-specified 
% scalar or a matching dimension numeric array. Syntax:
%
%                        C=times(A,B)
%
% Parameters:
%
%      A - a matrix or a cell array thereof
%
%      B - a matrix or a cell array thereof
%
% Outputs:
%
%      C - the resulting cell array
%
% Note: both arguments cannot be cell arrays.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cell/times.m>

function C=times(A,B)

% Check consistency
grumble(A,B)

% Decide the topology
if iscell(A)&&isnumeric(B)
    
    % Multiply every cell from the left
    for n=1:numel(A)
        A{n}=A{n}*B(n);
    end
    C=A;
    
elseif isnumeric(A)&&iscell(B)
    
    % Multiply every cell from the right
    for n=1:numel(B)
        B{n}=A(n)*B{n};
    end
    C=B;
    
else
    
    % Complain and bomb out
    error('at least one argument must be numeric.');
    
end

end

% Consistency enforcement
function grumble(A,B)
if (~iscell(A))&&(~isnumeric(A))
    error('A must be either numeric or a cell array.');
end
if (~iscell(B))&&(~isnumeric(B))
    error('B must be either numeric or a cell array.');
end
if ~all(size(A)==size(B),'all')
    error('the array of scalars and the cell array must have the same dimensions.');
end
end

% One of the great mistakes is to judge policies and programs 
% by their intentions rather than their results.
%
% Milton Friedman

