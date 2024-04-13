% Subtracts cell arrays element-by-element. Syntax:
%
%                       A=minus(A,B)
%
% Parameters:
%
%         A,B     - cell arrays of identical topology
%
% Outputs:
%
%         A       - the resulting cell array
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cell/minus.m>

function C=minus(A,B)

% Check consistency
grumble(A,B)

% Decide the topology
if iscell(A)&&iscell(B)
    
    % Subtract cell-by-cell
    for n=1:numel(A)
        A{n}=A{n}-B{n};
    end
    C=A;
    
elseif iscell(A)&&isnumeric(B)
    
    % Subtract from each cell
    for n=1:numel(A)
        A{n}=A{n}-B;
    end
    C=A;
    
elseif isnumeric(A)&&iscell(B)
    
    % Subtract from each cell
    for n=1:numel(B)
        B{n}=A-B{n};
    end
    C=B;
    
else
    
    % Complain and bomb out
    error('unsupported argument combination.');
    
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

% I can, therefore I am.
%
% Simone Weil

