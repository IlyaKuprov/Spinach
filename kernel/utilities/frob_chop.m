% Truncates SVD decomposition to the user-specified threshold
% in the Frobenius norm. Syntax:
%
%                     r=frob_chop(s,tol)
%
% Parameters:
%
%    s   - a vector of singular values for a matrix,
%          in descending order
%
%    tol - truncation threshold
%
% Outputs:
%
%    r   - the number of singular values to keep
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=frob_chop.m>

function r=frob_chop(s,tol)

% Check consistency
grumble(s,tol);

% Find the cutting point
x=cumsum(s(end:-1:1).^2);
k=find(x>=tol^2,1);

% Treat the zero case
if isempty(k)
    r=0;
else
    r=numel(s)-k+1;
end

end

% Consistency enforcement
function grumble(s,tol)
if (~isnumeric(tol))||(~isreal(tol))||(~isscalar(tol))||(tol<0)
    error('tol must be a non-negative real scalar.');
end
if (~isnumeric(s))||(~isreal(s))||(~isvector(s))||any(s<0)
    error('s must be a vector of non-negative real numbers.');
end
end

% "Morally equal": not equal, but must be treated as such.

