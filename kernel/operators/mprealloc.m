% Preallocates an operator in the current basis. Syntax:
%
%             A=mprealloc(spin_system,nnzpc)
%
% Parameters:
%
%    nnzpc  -  expected number of non-zeros per column
%
% Outputs:
%
%    A      -  all-zero sparse matrix of the appropriate
%              dimension with room for the specified num-
%              ber of non-zeroes
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=mprealloc.m>

function A=mprealloc(spin_system,nnzpc)

% Check consistency
grumble(spin_system,nnzpc);

% Preallocate the direct-sum operator
problem_dim=spin_system.bas.offsets(end);
A=spalloc(problem_dim,problem_dim,nnzpc*problem_dim);

end

% Consistency enforcement
function grumble(spin_system,nnzpc)
if (~isnumeric(nnzpc))||(~isreal(nnzpc))||(~isscalar(nnzpc))||(mod(nnzpc,1)~=0)
    error('nnzpc parameter must be a positive real integer.');
end
if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))
    error('the spin_system object does not contain the necessary data.');
end
end

% In 2006, Oxford's Magdalen College, where Erwin Schrodinger was a Fellow
% between 1933 and 1936, received a sum of money from a benefactor towards
% "increasing the art content of the College". A number of works were pre-
% sented for competition, among them a beautiful stone obelisk, called "Mo-
% nument to Knowledge", with Schrodinger's equation inscribed on it. The
% obelisk was rejected - the inscriber had missed the bar off Planck's con-
% stant. As the College Governing Body put it, the equation, as written,
% "would have exploded the stone it was inscribed upon".

