% Horizontal concatenation for RCV objects. Syntax:
%
%                     obj=horzcat(obj,obj2)
%
% Parameters:
%
%    obj   - left RCV object
%    obj2  - right RCV object
%
% Outputs:
%
%    obj   - concatenated object representing [obj obj2]
%
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rcv/horzcat.m>

function obj=horzcat(obj,obj2)

% Check consistency
grumble(obj,obj2);

% Align data locations between operands
if obj.isGPU&&~obj2.isGPU
    obj2.row=gpuArray(obj2.row);
    obj2.col=gpuArray(obj2.col);
    obj2.val=gpuArray(obj2.val);
    obj2.isGPU=true;
elseif ~obj.isGPU&&obj2.isGPU
    obj.row=gpuArray(obj.row);
    obj.col=gpuArray(obj.col);
    obj.val=gpuArray(obj.val);
    obj.isGPU=true;
end

% Shift column indices of the right operand
obj2.col=obj2.col+obj.numCols;

% Concatenate row, column and value arrays
obj.row=[obj.row;obj2.row];
obj.col=[obj.col;obj2.col];
obj.val=[obj.val;obj2.val];

% Update column count
obj.numCols=obj.numCols+obj2.numCols;

end

% Consistency enforcement
function grumble(obj,obj2)
if ~isa(obj,'rcv')||~isa(obj2,'rcv')
    error('both inputs must be rcv objects.');
end
if ~isscalar(obj)||~isscalar(obj2)
    error('concatenation inputs must be scalar rcv objects.');
end
if obj.numRows~=obj2.numRows
    error('row counts must match for horizontal concatenation.');
end
end

% The back half of your forties is a cursed age. It's not
% so much that (as the cliché goes) the policemen look
% younger - it's more that you're gloomily aware that
% you haven't done anything that they'd even consider
% arresting you for.
%
% Sam Leith
