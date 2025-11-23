% Returns the number of non-zero elements of the RCV object
%
% m.keitel@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function n = nnz(A)
    
    % Count unique index pairs in the
    % row and column index arrays
    n=numunique([A.row A.col],'rows');
    
end

% "By accepting insults and expressing gratitute for them."
%
% An old courtier, quoted by Seneca, 
% when asked how he had lasted so 
% long in the imperial service.

