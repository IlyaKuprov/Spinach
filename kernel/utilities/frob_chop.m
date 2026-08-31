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
% Note: singular values below numel(s)*eps*max(s), the standard
%       numerical rank threshold, are set to zero because the
%       SVD that produced them does not resolve them; above that
%       threshold the requested tolerance is honoured exactly.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=frob_chop.m>

function r=frob_chop(s,tol)

% Check consistency
grumble(s,tol);

% Remove round-off artefacts below the numerical rank threshold
s=real(s(:));
s(abs(s)<numel(s)*eps*max(abs(s)))=0;
s=max(s,0);

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
if (~isnumeric(s))||((~isvector(s))&&(~isempty(s)))
    error('s must be a vector of non-negative real numbers.');
end
if any(~isfinite(s(:)))||any(abs(imag(s(:)))>1e-10*max(abs(s(:))))||...
   any(real(s(:))<-1e-10*max(abs(s(:))))
    error('s must be a vector of non-negative real numbers.');
end
end

% "Morally equal": not equal, but must be treated as such.

