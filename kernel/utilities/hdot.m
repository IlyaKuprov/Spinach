% Hadamard matrix product. Useful as a replacement for 
% trace(A'*B) because trace(A'*B)=hadm(conj(A),B) and 
% the latter only needs O(n^2) multiplications as oppo-
% sed to O(n^3) for trace(A'*B). Syntax: 
%
%                     H=hdot(A,B)
%
% where A and B are square matrices and H is a scalar.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hdot.m>

function H=hdot(A,B)

% Check consistency
grumble(A,B);

% Do the calculation
H=sum(conj(A).*B,'all');

end

% Consistency enforcement
function grumble(A,B)
if (~isnumeric(A))||(~isnumeric(B))
    error('both inputs must be numeric.');
end
if ~all(size(A)==size(B))
    error('the two inputs must have identical dimensions.');
end
end

% An infinite number of mathematicians walk into a bar. The first one
% orders a pint of beer, the second one half a pint, the third one a
% quarter... "Gotcha!" says the bartender and pours two pints.

