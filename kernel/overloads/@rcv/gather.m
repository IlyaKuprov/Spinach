% Gathers an RCV object from the GPU back to the CPU. Syntax:
%
%                       obj=gather(obj)
%
% Parameters:
%
%    obj   - an RCV object
%
% Outputs:
%
%    obj   - the same object with data stored on the CPU
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/gather.m>

function obj=gather(obj)

% Check consistency
grumble(obj);

% Move row, column and value data to the CPU
if obj.isGPU
    obj.row=gather(obj.row);
    obj.col=gather(obj.col);
    obj.val=gather(obj.val);
    obj.isGPU=false;
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

% Aerie, I've noticed the unfortunate fact that you live
% by one of the great lessons of history that nothing is
% often a good thing to do and a clever thing to say.
%
% Edwin Odesseiron, in Baldur's Gate 2
