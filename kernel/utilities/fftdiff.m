% Spectral differentiation kernel. Syntax:
%
%               kern=fftdiff(order,npoints,dx)
%
% Parameters:
%
%       order   - order of the derivative
%
%       npoints - number of points in the grid
%
%       dx      - grid step length
%
% Output:
%
%       kern - a vector that is to be used for accurate
%              numerical differentiation of periodic real
%              signals in the following way:
%
%               derivative=real(ifft(fft(signal).*kern));
%
% Note: periodic boundary conditions.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fftdiff.m>

function kern=fftdiff(order,npoints,dx)

% Check consistency
grumble(order,npoints,dx);

% Adapt to the point count
if mod(npoints,2)==1
    
    % Kernel for odd point counts
    kern=ifftshift((2i*pi*((1-npoints)/2:((npoints)/2))/(npoints*dx)).^order);
    
elseif mod(npoints,2)==0
    
    % Kernel for even point counts
    kern=ifftshift((2i*pi*(((-npoints)/2):((npoints-1)/2))/(npoints*dx)).^order);
    
else
    
    % Complain and bomb out
    error('npoints parameter must be an integer.');
    
end

end

% Consistency enforcement
function grumble(order,npoints,dx)
if (~isnumeric(order))||(~isreal(order))||(numel(order)~=1)||...
   (order<1)||(mod(order,1)~=0)
    error('order must be a non-negative real integer.');
end
if (~isnumeric(npoints))||(~isreal(npoints))||(numel(npoints)~=1)||...
   (npoints<1)||(mod(npoints,1)~=0)
    error('npoints must be a non-negative real integer.');
end
if (~isnumeric(dx))||(~isreal(dx))||(numel(dx)~=1)||(dx<=0)
    error('dx must be a positive real number.');
end
end

% A few people laughed, a few people cried, most people were silent. I re-
% membered the line from the Hindu scripture, the Bhagavad-Gita; Vishnu is
% trying to persuade the Prince that he should do his duty and, to impress
% him, takes on his multi-armed form and says, 'Now I am become Death, the
% destroyer of worlds.' I suppose we all thought that, one way or another.
%
% J. Robert Oppenheimer, about the first atomic detonation

