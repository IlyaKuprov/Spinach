% RCV object size.
%
% m.keitel@soton.ac.uk

function s=size(obj,dim)
    if nargin==1
        s=[obj.numRows obj.numCols]; % Return both dimensions
    else
        if dim==1
            s=obj.numRows; % Return number of rows
        elseif dim==2
            s=obj.numCols; % Return number of columns
        else
            s=1; % Consistent with MATLAB behavior for higher dimensions
        end
    end
end

% I did not succeed in life by intelligence. I succeeded
% because I have a long attention span.
%
% Charlie Munger

