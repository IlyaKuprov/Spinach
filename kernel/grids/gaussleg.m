% Computes Gauss-Legendre points and weights in [a,b] interval
% with accuracy order n. Syntax:
%
%                     [x,w]=gaussleg(a,b,n)
%
% Parameters:
%
%    a - left edge of the interval
%
%    b - right edge of the interval
%
%    n - accuracy order, the number of points in the 
%        resulting grid will be n+1.
%
% Outputs:
%
%    x - Gauss-Legendre points
%
%    w - Gauss-Legendre weights
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=gaussleg.m>

function [x,w]=gaussleg(a,b,n)

% Check consistency
grumble(a,b,n);

% Initial guess for the nodes in [-1 1]
x=cos((2*(0:n)'+1)*pi/(2*n+2))+(0.27/(n+1))*sin(pi*linspace(-1,1,n+1)'*n/(n+2));

% Newton-Raphson refinement
V=zeros(n+1,n+2); prev_x=0;
while max(abs(x-prev_x))>eps
    V(:,1)=1; V(:,2)=x;
    for k=2:(n+1)
        V(:,k+1)=((2*k-1)*x.*V(:,k)-(k-1)*V(:,k-1))/k;
    end
    dV=(n+2)*(V(:,n+1)-x.*V(:,n+2))./(1-x.^2);
    prev_x=x; x=x-V(:,n+2)./dV;
end

% Compute the weights
w=(b-a)./((1-x.^2).*dV.^2)*((n+2)/(n+1))^2;

% Map from [-1,1] to [a,b]
x=(a*(1-x)+b*(1+x))/2;

end

% Consistency enforcement
function grumble(a,b,n)
if (~isnumeric(a))||(~isreal(a))||(~isfinite(a))||(numel(a)~=1)
    error('a must be a finite real number.');
end
if (~isnumeric(b))||(~isreal(b))||(~isfinite(b))||(numel(b)~=1)
    error('b must be a finite real number.');
end
if a >= b
    error('a must be smaller than b.');
end
if (~isnumeric(n))||(~isreal(n))||(~isfinite(n))||...
   (numel(b)~=1)||(n<1)||mod(n,1)
    error('n must be a positive real integer.');
end
if n>50
    error('GL schemes with n>50 are unstable - subdivide your interval.');
end
end

% Life is a fountain of delight; but where the rabble also 
% drinks all wells are poisoned.
%
% Friedrich Nietzsche

