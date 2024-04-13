% Subtracts an appropriate multiple of the unit matrix to make 
% the input matrix traceless. Syntax:
%
%                        A=remtrace(A)
%
% Parameters:
%
%   A - a square matrix
%
% Outputs:
%
%   A - a square matrix with a zero trace
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=remtrace.m>

function A=remtrace(A)

% Check consistency
grumble(A);

% Kill the trace
dim=size(A,1); A=A-speye(dim)*trace(A)/dim;

end

% Consistency enforcement
function grumble(A)
if (~isnumeric(A))||(size(A,1)~=size(A,2))
    error('A must be a square matrix.');
end
end

% Cacophobia, n.
%
% The fear of ugliness and of things
% that are ugly.

