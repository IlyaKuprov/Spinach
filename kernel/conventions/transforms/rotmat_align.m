% Rotation matrix aligning two vectors.
%
% alexey.bogdanov@weizmann.ac.il

function rot_mat=rotmat_align(v_from,v_to)

% Normalize the input vectors
v_from=v_from/norm(v_from);
v_to=v_to/norm(v_to);

% Compute the rotation axis and angle
axis_vec=cross(v_from,v_to);
sin_ang=norm(axis_vec);
cos_ang=dot(v_from,v_to);

% Handle parallel vectors
if sin_ang<1e-12
    if cos_ang>0
        rot_mat=eye(3);
        return
    else
        axis_vec=null(v_from.');
        axis_vec=axis_vec(:,1);
        sin_ang=0;
        cos_ang=-1;
    end
end

% Normalize the rotation axis
axis_vec=axis_vec/sin_ang;

% Build the cross-product matrix
skew=[0           -axis_vec(3)  axis_vec(2);
      axis_vec(3)  0           -axis_vec(1);
     -axis_vec(2)  axis_vec(1)  0         ];

% Build the rotation matrix
rot_mat=eye(3)+sin_ang*skew+(1-cos_ang)*(skew*skew);

end

% Dignity does not consist in possessing 
% honours, but in deserving them.
%
% Aristotle

