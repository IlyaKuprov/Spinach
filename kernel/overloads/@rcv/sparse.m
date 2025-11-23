% Convert to sparse matrix.
%
% m.keitel@soton.ac.uk

function S=sparse(obj)

    if isempty(obj.col)

        % Empty sparse matrix
        S=spalloc(obj.numRows,obj.numCols,0);

    else

        % Matlab's constructor
        S=sparse(obj.row,obj.col,obj.val,obj.numRows,obj.numCols);

    end
    
end

% Working 16 hours a day, 7 days a week, 52 weeks
% in a year, and people still calling me lucky.
%
% Elon Musk

