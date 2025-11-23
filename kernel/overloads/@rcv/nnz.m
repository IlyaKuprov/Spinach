% Returns the number of non-zero elements in an RCV object. Syntax:
%
%                          n=nnz(A)
%
% Parameters:
%
%    A     - RCV object
%
% Outputs:
%
%    n     - number of non-zero entries
%
% m.keitel@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rcv/nnz.m>

function n=nnz(A)

% Check consistency
grumble(A);

% Count unique index pairs
n=numunique([A.row A.col],'rows');

end

% Consistency enforcement
function grumble(A)
if ~isa(A,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% "By accepting insults and expressing gratitute for them."
%
% An old courtier, quoted by Seneca,
% when asked how he had lasted so
% long in the imperial service.

