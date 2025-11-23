% Transfers object ot the GPU
%
% m.keitel@soton.ac.uk

function obj=gpuArray(obj)

if ~obj.isGPU
    obj.row=gpuArray(obj.row);
    obj.col=gpuArray(obj.col);
    obj.val=gpuArray(obj.val);
    obj.isGPU=true;
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

