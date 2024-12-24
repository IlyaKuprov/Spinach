% Analytical Tikhonov regularised solution to K*x=y without any
% constraints (indeterminate output). Syntax:
%
%                [x,err,reg]=tikhoind(K,D,y,lam)
%
% Parameters:
%
%    K      - kernel matrix, may be complex, may be non-square
%
%    D      - regularisation matrix
%
%    y      - a column vector, may be complex
%
%    lam    - Tikhonov regularisation parameter
%
% Outputs:
%
%    x      - a real vector, a minimum of 
%             norm(K*x-y,2)^2+lambda*norm(D*x,2)^2
%
%    err    - error signal norm(K*x-y,2)^2
%
%    reg    - regularisation signal norm(D*x,2)^2
%
% Note: for best numerical performance, scale K to have approxima-
%       tely unit 2-norm, and y to have approximately unit 1-norm.
%
% Note: see tikhonov.m for the positive-constraned solver.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=tikhoind.m>

function [x,err,reg]=tikhoind(K,D,y,lam)

% Check consistency
grumble(K,D,y,lam);

% Analytical solution
x=((K'*K)+lam*(D'*D))\(K'*y);

% Error and regularisation signals
if nargout>1, err=norm(K*x-y,2)^2; end
if nargout>2, reg=norm(D*x,2)^2; end

end

% Consistency enforcement
function grumble(K,D,y,lam)
if (~isnumeric(K))||(~isnumeric(D))||...
   (~isnumeric(y))||(~isnumeric(lam))
    error('all inputs must be numeric.');
end
if size(K,1)~=size(y,1)
    error('dimensions of K and y are not consistent.');
end
if (~isreal(lam))||(~isscalar(lam))||(lam<0)
    error('lam must be a positive real scalar.');
end
end

% It is easy to love a country that is famous for 
% chocolate and beer.
% 
% Barack Obama, about Belgium

