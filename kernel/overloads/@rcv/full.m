% Convert to full matrix
%
% m.keitel@soton.ac.uk

function S=full(obj)

    % Generate the sparse matrix and then make it full
    S=full(sparse(obj));

end

% Whenever you find yourself on the side of the 
% majority, it is time to pause and reflect.
%
% Mark Twain

