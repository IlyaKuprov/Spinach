% Non-conjugating transpose of a matrix-free core. Syntax:
%
%                         core=transpose(core)
%
% Parameters:
%
%    core - a matfree core
%
% Outputs:
%
%    core - transposed core with swapped dimensions
%
% ilya.kuprov@weizmann.ac.il

function core=transpose(core)

% Check consistency
grumble(core);

% Transpose through conjugated forward and adjoint actions
forward=core.forward; adjoint=core.adjoint;
core.forward=@(block)conj(adjoint(conj(block)));
core.adjoint=@(block)conj(forward(conj(block)));
core.dims=fliplr(core.dims);

end

% Consistency enforcement
function grumble(core)
if ~isa(core,'matfree'), error('core must be a matfree object.'); end
end


