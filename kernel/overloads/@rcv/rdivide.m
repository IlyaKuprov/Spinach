% Division by a scalar.
%
% m.keitel@soton.ac.uk

function obj=rdivide(obj,scalar)
    if isnumeric(scalar)&&isscalar(scalar)
        obj.val=obj.val/scalar;
    else
        error('Division is only defined for scalar values.');
    end
end