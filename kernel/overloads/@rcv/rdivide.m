% Divides an RCV object by a numeric scalar. Syntax:
%
%                      obj=rdivide(obj,scalar)
%
% Parameters:
%
%    obj    - RCV object
%    scalar - numeric scalar
%
% Outputs:
%
%    obj    - result of obj./scalar
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/rdivide.m>

function obj=rdivide(obj,scalar)

% Check consistency
grumble(obj,scalar);

% Divide stored values by the scalar
obj.val=obj.val/scalar;

end

% Consistency enforcement
function grumble(obj,scalar)
if ~isa(obj,'rcv')
    error('the first argument must be an rcv object.');
end
if ~isscalar(obj)
    error('the rcv input must be a scalar object.');
end
if ~(isnumeric(scalar)&&isscalar(scalar))
    error('division is only defined for numeric scalars.');
end
end
