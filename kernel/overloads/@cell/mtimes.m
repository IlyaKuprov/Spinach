% Multiplies all entries of a cell array by a user-specified 
% matrix. Syntax:
%
%                        C=mtimes(A,B)
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
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cell/mtimes.m>

function C=mtimes(A,B)

% Check consistency
grumble(A,B)

% Decide the topology
if iscell(A)&&isnumeric(B)
    
    % Multiply every cell from the left
    for n=1:numel(A)
        A{n}=A{n}*B;
    end
    C=A;
    
elseif isnumeric(A)&&iscell(B)
    
    % Multiply every cell from the right
    for n=1:numel(B)
        B{n}=A*B{n};
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
end

% We keep blaming Comrade Stalin and I agree that we have 
% good reasons. But I would like to ask: who wrote those
% four million denunciation letters?
% 
% S.D. Dovlatov

