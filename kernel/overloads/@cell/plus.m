% Adds cell arrays element-by-element. Syntax:
%
%                       A=plus(A,B)
%
% Parameters:
%
%         A,B     - cell arrays of identical topology
%
% Outputs:
%
%         C       - the resulting cell array
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cell/plus.m>

function C=plus(A,B)

% Check consistency
grumble(A,B)

% Decide the topology
if iscell(A)&&iscell(B)
    
    % Add cell-by-cell
    for n=1:numel(A)
        A{n}=A{n}+B{n};
    end
    C=A;
    
elseif iscell(A)&&isnumeric(B)
    
    % Add to each cell
    for n=1:numel(A)
        A{n}=A{n}+B;
    end
    C=A;
    
elseif isnumeric(A)&&iscell(B)
    
    % Add to each cell
    for n=1:numel(B)
        B{n}=A+B{n};
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

% I came into the room, which was half dark, and presently spotted Lord
% Kelvin in the audience and realized that I was in for trouble at the last
% part of my speech dealing with the age of the Earth, where my views
% conflicted with his. To my relief, Kelvin fell fast asleep, but as I
% came to the important point, I saw the old bird sit up, open an eye and
% cock a baleful glance at me! Then a sudden inspiration came, and I said
% "Lord Kelvin had limited the age of the Earth, provided no new source of
% energy was discovered. That prophetic utterance refers to what we are
% now considering tonight, radium." Behold! The old boy beamed upon me.
% 
% Ernest Rutherford

