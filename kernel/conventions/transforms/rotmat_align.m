% Rotation matrix aligning one vector with another vector. Syntax:
%
%                   rot_mat=rotmat_align(v_from,v_to)
%
% Parameters:
%
%    v_from   - three-element real vector to rotate
%
%    v_to     - three-element real vector to align to
%
% Outputs:
%
%    rot_mat  - 3x3 rotation matrix that satisfies
%               rot_mat*(v_from/norm(v_from,2))=v_to/norm(v_to,2)
%
% Note: aligning one vector with another leaves one rotational degree of
%       freedom around the aligned direction. If the resulting matrix is
%       converted to ZYZ Euler angles, this freedom appears as a non-uni-
%       que third Euler angle (the twist around the aligned axis). This
%       implementation fixes that freedom by returning the minimum-
%       angle alignment without any additional twist around the aligned
%       direction. In the anti-parallel case the axis is non-unique, and
%       the first null-space basis vector orthogonal to v_from is used.
%
% alexey.bogdanov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rotmat_align.m>

function rot_mat=rotmat_align(v_from,v_to)

% Check consistency
grumble(v_from,v_to);

% Normalise the input vectors
v_from=v_from(:)/norm(v_from,2);
v_to=v_to(:)/norm(v_to,2);

% Set the collinearity tolerance
align_tol=1e-12;

% Compute the alignment axis and the alignment angle cosine
rot_axis=cross(v_from,v_to);
axis_norm=norm(rot_axis,2);
cos_ang=dot(v_from,v_to);

% Handle parallel and anti-parallel vectors
if axis_norm<align_tol
    if cos_ang>0

        % Return identity when no rotation is needed
        rot_mat=eye(3);
        return

    else

        % Choose an arbitrary axis orthogonal to the input vector
        null_basis=null(v_from.');
        rot_axis=null_basis(:,1);
        axis_norm=norm(rot_axis,2);

        % Set the anti-parallel rotation angle parameters
        sin_ang=0; cos_ang=-1;

    end
else

    % Use the axis norm as the sine of the rotation angle
    sin_ang=axis_norm;

end

% Normalise the rotation axis
rot_axis=rot_axis/axis_norm;

% Build the cross-product matrix of the rotation axis
skew_mat=[0            -rot_axis(3)  rot_axis(2);
          rot_axis(3)   0           -rot_axis(1);
         -rot_axis(2)   rot_axis(1)  0          ];

% Build the rotation matrix using Rodrigues formula
rot_mat=eye(3)+sin_ang*skew_mat+(1-cos_ang)*(skew_mat*skew_mat);

end

% Consistency enforcement
function grumble(v_from,v_to)
if (~isnumeric(v_from))||(~isreal(v_from))||(~isvector(v_from))||...
   (numel(v_from)~=3)||any(~isfinite(v_from(:)))
    error('v_from must be a real finite vector with three elements.');
end
if (~isnumeric(v_to))||(~isreal(v_to))||(~isvector(v_to))||...
   (numel(v_to)~=3)||any(~isfinite(v_to(:)))
    error('v_to must be a real finite vector with three elements.');
end
if norm(v_from,2)<eps('double')
    error('v_from must not be a zero vector.');
end
if norm(v_to,2)<eps('double')
    error('v_to must not be a zero vector.');
end
end

% Dignity does not consist in possessing 
% honours, but in deserving them.
%
% Aristotle


