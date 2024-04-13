% Normalized Gaussian function in magnetic resonance 
% notation. Syntax:
%
%                g=gaussfun(x,fwhm)
%
% Parameters:
%
%        x - argument values, a real array of any dimension
%
%     fwhm - full width at half-maximum
%
% Outputs:
%
%        g - function values at the points specified in x
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=gaussfun.m>

function y=gaussfun(x,fwhm)

% Check consistency
grumble(x,fwhm);

% Compute standard deviation
sigma=fwhm/(2*sqrt(2*log(2)));

% Compute the Gaussian
y=(1/(sigma*sqrt(2*pi)))*exp(-(x.^2)/(2*sigma^2));

end

% Consistency enforcement
function grumble(x,fwhm)
if (~isnumeric(x))||(~isreal(x))
    error('x must be an array of real numbers.');
end
if (~isnumeric(fwhm))||(~isreal(fwhm))||...
   (numel(fwhm)~=1)||(fwhm<=0)
    error('fwhm must be a positive real number.');
end
end

% Fifty years ago the back streets of Leningrad
% have taught me one lesson: when a fight is un-
% avoidable, punch first.
%
% Vladimir Putin

