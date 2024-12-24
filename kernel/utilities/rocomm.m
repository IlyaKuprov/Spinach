% Right-ordered nested commutator [[[[A{1},A{2}],A{3}],A{4}],...] 
% built from the user-supplied matrices. Syntax:
%
%                           C=rocomm(A)
%
% Parameters:
%
%     A   - a cell array of square matrices
%
% Outputs:
%
%     C   - right-ordered nested commutator
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rocomm.m>

function C=rocomm(A)

% Check consistency
grumble(A);

% Nest the commutators
C=A{1};
for n=2:numel(A)
    C=C*A{n}-A{n}*C;
end

end

% Consistency enforcement
function grumble(A)
if (~iscell(A))||(~all(cellfun(@isnumeric,A),'all'))||...
   (~all(cellfun(@(x)size(x,1),A)==cellfun(@(x)size(x,2),A),'all'))
    error('A must be a cell array of square matrices.');
end
end

% It's not science I don't trust - it's the scientists.
%
% James Delingpole, 
% a climate sceptic

