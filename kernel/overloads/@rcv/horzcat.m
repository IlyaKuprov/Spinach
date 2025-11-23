% Horizontal Concatenation [A, B]
%
% m.keitel@soton.ac.uk

function obj=horzcat(obj,obj2)
    if obj.numRows~=obj2.numRows
        error('Row counts must match for horizontal concatenation.');
    end
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
    obj2.col=obj2.col+obj.numCols;
    obj.row=[obj.row; obj2.row];
    obj.col=[obj.col; obj2.col];
    obj.val=[obj.val; obj2.val];
    obj.numCols=obj.numCols+obj2.numCols;
end

% The back half of your forties is a cursed age. It's not 
% so much that (as the cliché goes) the policemen look 
% younger - it's more that you're gloomily aware that
% you haven't done anything that they'd even consider
% arresting you for.
%
% Sam Leith

