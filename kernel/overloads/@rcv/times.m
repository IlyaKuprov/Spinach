% Multiplies an RCV sparse matrix by a numeric scalar. Syntax:
%
%                      obj=times(scalar,obj)
%
% Parameters:
%
%    scalar - numeric scalar
%
%    obj    - RCV sparse matrix
%
% Outputs:
%
%    obj    - RCV sparse matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/times.m>

function obj=times(scalar,obj)

% Check consistency
grumble(scalar,obj);

% Multiply values by the scalar
obj.val=obj.val*scalar;

end

% Consistency enforcement
function grumble(scalar,obj)
if ~isa(obj,'rcv')
    error('the second argument must be an RCV sparse matrix.');
end
if (~isnumeric(scalar))||(~isscalar(scalar))
    error('multiplication is only defined by numeric scalars.');
end
end

% They say that the fish that gets away
% looks bigger than it really is.
%
% Seven Samurai film (1954)

