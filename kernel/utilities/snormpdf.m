% Azzalini's skew normal distribution. Syntax:
%
%              p=snormpdf(x,mu,sigma,alpha)
%
% Parameters:
%
%   x     - an array of real numbers
%
%   mu    - expectation value of the normal distribution
%
%   sigma - standard deviation of the normal distribution
%
%   alpha - skew factor, a real number
%
% Outputs:
%
%   p     - an array of probability densities, 
%           same shape as x 
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=snormpdf.m>

function p=snormpdf(x,mu,sigma,alpha)

% Check consistency
grumble(x,mu,sigma,alpha);

% Equation 2 in http://www.jstor.org/stable/4615982
p=2*normpdf(x,mu,sigma).*normcdf(alpha*x,alpha*mu,sigma);

end

% Consistency enforcement
function grumble(x,mu,sigma,alpha)
if (~isnumeric(x))||(~isreal(x))
    error('x must be a real numeric array.');
end
if (~isnumeric(mu))||(~isreal(mu))||(~isscalar(mu))
    error('mu must be a real scalar.');
end
if (~isnumeric(sigma))||(~isreal(sigma))||...
   (~isscalar(sigma))||(sigma<=0)
    error('sigma must be a real positive scalar.');
end
if (~isnumeric(alpha))||(~isreal(alpha))||(~isscalar(alpha))
    error('alpha must be a real scalar.');
end
end

% The smallest minority on earth is the individual. Those who deny
% individual rights cannot claim to be defenders of minorities.
%
% Ayn Rand

