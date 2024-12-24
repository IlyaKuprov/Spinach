% Logarithm of the factorial function. Avoids complications with
% factorials of large numbers overflowing 64-bit numbers. Syntax:
%
%                       lf=logfactorial(n)
%
% Parameters:
%
%        n - non-negative integer number
%
% Outputs:
%
%       lf - logarithm of the factorial of n
%
% Notes: double precision overflow is a persistent problem with
%        Clebsch-Gordan coefficients and other objects that in-
%        volve factorials.
% 
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=logfactorial.m>

function lf=logfactorial(n)

% Check consistency
grumble(n);

% Use the built-in log(gamma(n)) function
lf=gammaln(n+1);

end

% Consistency enforcement
function grumble(n)
if (~isnumeric(n))||(~isreal(n))||...
   any(n(:)<0)||any(mod(n(:),1)~=0)
    error('elements of n must be non-negative integers.');
end
end

% For ten years or so, my name was "that jerk". But 
% that was a promotion. Before, I was "Who's he?"
%
% Richard Thaler

