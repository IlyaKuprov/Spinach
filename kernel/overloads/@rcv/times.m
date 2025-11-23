% Multiplication by a scalar.
%
% m.keitel@soton.ac.uk

function obj=times(scalar,obj)
    if isnumeric(scalar)&&isscalar(scalar)
        obj.val=obj.val*scalar;
    else
        error('Multiplication is only defined for scalar values.');
    end
end   

% They say that the fish that gets away
% looks bigger than it really is.
%
% Seven Samurai film (1954)

