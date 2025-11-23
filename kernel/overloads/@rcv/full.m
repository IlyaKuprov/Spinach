% Converts an RCV sparse matrix into a full matrix. Syntax:
%
%                       S=full(obj)
%
% Parameters:
%
%    obj   - an RCV sparse matrix
%
% Outputs:
%
%    S     - a full matrix
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/full.m>

function S=full(obj)

% Check consistency
grumble(obj);

% Delegate to Matlab
S=full(sparse(obj));

end

% Consistency enforcement
function grumble(obj)
if ~isa(obj,'rcv')
    error('the input must be an RCV sparse matrix.');
end
end

% Whenever you find yourself on the side of the
% majority, it is time to pause and reflect.
%
% Mark Twain

