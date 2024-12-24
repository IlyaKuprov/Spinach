% Exponential drop function. Produces an exponential fall-off
% from a specified value to a specified value with the specified
% rate and the number of points. Syntax:
%
%       drop=expdrop(from,to,duration,npoints,drop_rate)
%
% Parameters:
%
%        from - the value to drop from
%
%          to - the value to drop to
%
%    duration - drop duration, seconds
%
%     npoints - the number of discretisation points
%               in the drop
%
%   drop_rate - exponential drop rate, Hz
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=expdrop.m>

function drop=expdrop(from,to,duration,npoints,drop_rate)

% Check consistency
grumble(from,to,duration,npoints,drop_rate);

% Get the exponential drop parameters
B=(from-to)/(1-exp(-drop_rate*duration)); A=from-B;

% Compute the drop
drop=A+B*exp(-drop_rate*linspace(0,duration,npoints));

end

% Consistency enforcement
function grumble(from,to,duration,npoints,drop_rate)
if (~isnumeric(npoints))||(~isreal(npoints))||...
   (numel(npoints)~=1)||(npoints<1)||(mod(npoints,1)~=0)
    error('npoints variable must be a positive real integer');
end
if (~isnumeric(from))||(~isreal(from))||(numel(from)~=1)
    error('from variable must be a real number');
end
if (~isnumeric(to))||(~isreal(to))||(numel(to)~=1)
    error('to variable must be a real number');
end
if (~isnumeric(duration))||(~isreal(duration))||...
   (numel(duration)~=1)||(duration<=0)
    error('duration variable must be a positive real number');
end
if (~isnumeric(drop_rate))||(~isreal(drop_rate))||...
   (numel(drop_rate)~=1)||(drop_rate<=0)
    error('drop_rate variable must be a positive real number');
end
end

% If you could reason with religious people, there
% would be no religious people.
%
% Gregory House

