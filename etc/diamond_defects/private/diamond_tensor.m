% Tensor builder for diamond-defect principal values. Syntax:
%
%          M=diamond_tensor(values,frame)
%
% This internal helper builds a symmetric tensor from principal values
% and a principal-axis frame.
%
% Parameters:
%
%    values  - three principal values
%
%    frame   - three principal axes as columns
%
% Outputs:
%
%    M       - symmetric tensor matrix
%
% <https://spindynamics.org/wiki/index.php?title=diamond_tensor.m>

function M=diamond_tensor(values,frame)

% Check consistency
grumble(values,frame);

% Orthogonalise the principal-axis frame
frame=diamond_frame_orth(frame);

% Build and symmetrise the tensor
M=frame*diag(values)*frame';
M=(M+M')/2;

end

% Consistency enforcement
function grumble(values,frame)
if(~isnumeric(values)||~isreal(values)||numel(values)~=3)
    error('values must be a real three-element array.');
end
if(~isnumeric(frame)||~isreal(frame)||~isequal(size(frame),[3 3]))
    error('frame must be a real 3x3 matrix.');
end
end

% Tensor numbers are meaningless without their axes.

