% Hermitian conjugate of a matrix-free core. Syntax:
%
%                         core=ctranspose(core)
%
% Parameters:
%
%    core - a matfree core
%
% Outputs:
%
%    core - adjoint core with swapped dimensions and actions
%
% ilya.kuprov@weizmann.ac.il

function core=ctranspose(core)

% Check consistency
grumble(core);

% Swap the actions and conjugate the coefficient
forward=core.forward; core.forward=core.adjoint; core.adjoint=forward;
core.dims=fliplr(core.dims); core.coeff=conj(core.coeff);

end

% Consistency enforcement
function grumble(core)
if ~isa(core,'matfree'), error('core must be a matfree object.'); end
end


