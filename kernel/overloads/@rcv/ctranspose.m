% Returns conjugate transpose of RCV object
%
% m.keitel@soton.ac.uk

function obj=ctranspose(obj)

    % Efficiently swap rows and columns
    [obj.col,obj.row]=deal(obj.row,obj.col);

    % Conjugate values
    obj.val=conj(obj.val);

end

% At the Tower of London we remembered in 
% our prayers the elephant kept there by 
% James I, which, poor creature, was never
% given anything to drink but wine.
%
% Tom Holland

