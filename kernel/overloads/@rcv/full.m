% Converts an RCV object to a full matrix. Syntax:
%
%                       S=full(obj)
%
% Parameters:
%
%    obj   - an RCV object
%
% Outputs:
%
%    S     - a full matrix containing the data from obj
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/full.m>

function S=full(obj)

% Check consistency
grumble(obj);

% Convert to a sparse matrix and then expand to full
S=full(sparse(obj));

end

% Consistency enforcement
function grumble(obj)
if ~isa(obj,'rcv')
    error('the input must be an rcv object.');
end
if ~isscalar(obj)
    error('the input must be a scalar rcv object.');
end
end

% Whenever you find yourself on the side of the
% majority, it is time to pause and reflect.
%
% Mark Twain
