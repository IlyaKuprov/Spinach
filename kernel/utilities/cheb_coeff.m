% Discrete cosine transform algorithm for Chebyshev expansion
% coefficients of the user-specified scalar function. Syntax:
%
%                     c=cheb_coeff(f,a,b,n)
%
% Parameters:
%
%    f - function handle, must be vectorised
%
%    a - left edge of the expansion interval
%
%    b - right edge of the expansion interval
%
%    n - number of Chebyshev polynomials in 
%        the expansion
%
% Outputs:
%  
%    c - a vector of expansion coefficients
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cheb_coeff.m>

function c=cheb_coeff(f,a,b,n)

% Check consistency
grumble(f,a,b,n);

% [-1,+1] query points
x=cos(((1:n)*2-1)*pi/(2*n));

% Scaled query points
x=0.5*(a+x*(b-a)+b);

% Expansion coefficients
c=dct(f(x))/sqrt(n);
c(2:n)=c(2:n)*sqrt(2);

end

% Consistency enforcement
function grumble(f,a,b,n)
if ~isa(f,'function_handle')
    error('f must be a function handle.');
end
if (~isnumeric(a))||(~isscalar(a))||(~isreal(a))||...
   (~isnumeric(b))||(~isscalar(b))||(~isreal(b))||(a>=b)
    error('a and b must be real scalars such that a < b');
end
if (~isnumeric(n))||(~isscalar(n))||...
   (~isreal(n))||(mod(n,1)~=0)||(n<1)
    error('n must be a positive real integer.');
end
end

% Overheard in the Old Parsonage:
%
% Person A, horrified: "You would *do* that?!"
% Person B: "As a human, no. As a chemist..." [shrugs]

