% Transfers an RCV object to the GPU. Syntax:
%
%                     obj=gpuArray(obj)
%
% Parameters:
%
%    obj   - an RCV object
%
% Outputs:
%
%    obj   - the same object with data stored on the GPU
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/gpuArray.m>

function obj=gpuArray(obj)

% Check consistency
grumble(obj);

% Move row, column and value data to the GPU
if ~obj.isGPU
    obj.row=gpuArray(obj.row);
    obj.col=gpuArray(obj.col);
    obj.val=gpuArray(obj.val);
    obj.isGPU=true;
end

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

% Then it got worse. The book is very, very good. If
% someone's going to beat you to the punch with a great
% book idea, the least they can do is write something
% crap. Not Andrew. Which shouldn't really come as a
% surprise, since the little bastard is prodigiously
% talented.
%
% Toby Young, about Andrew Doyle's book titled
% "The New Puritans: How the Religion of Social
% Justice Captured the Western World"
