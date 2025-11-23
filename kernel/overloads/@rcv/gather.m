% Gathers an RCV sparse matrix from GPU. Syntax:
%
%                 obj=gather(obj)
%
% Parameters:
%
%    obj   - an RCV sparse matrix
%
% Outputs:
%
%    obj   - the same matrix with data stored on the CPU
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/gather.m>

function obj=gather(obj)

% Check consistency
grumble(obj);

% Gather to CPU
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
    error('the input must be an RCV sparse matrix.');
end
end

% Aerie, I've noticed the unfortunate fact that you live
% by one of the great lessons of history that nothing is
% often a good thing to do and a clever thing to say.
%
% Edwin Odesseiron, in Baldur's Gate 2

