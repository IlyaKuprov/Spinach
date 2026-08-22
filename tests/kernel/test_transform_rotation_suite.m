% Tests rotation transform helpers. Syntax:
%
%                    result=test_transform_rotation_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks active ZYZ Euler rotations, Euler/DCM inversion, rotation
% composition, angle-axis normalisation, quaternion round-trips, Wigner
% identity rotation, and minimum-angle vector alignment.
%
% ilya.kuprov@weizmann.ac.il

function result=test_transform_rotation_suite()

% Announce the test target
fprintf('TESTING: Rotation transform helpers\n');

% State the geometric target of the test
result=new_test_result('kernel/transform_rotation_suite',...
                       'Rotation transform helpers',...
                       'rotation transforms must preserve active rotation, composition, and round-trip conventions.');

% Check the active ZYZ convention on a known quarter-turn
R=euler2dcm(pi/2,0,0);
R_ref=[0 -1 0;1 0 0;0 0 1];
result=test_close(result,'euler2dcm active Z rotation',R,R_ref,1e-15,1e-15,...
                  'a positive active Z rotation maps the x axis into y');
result=test_close(result,'euler2dcm vector input',euler2dcm([pi/2 0 0]),R_ref,1e-15,1e-15,...
                  'the one-vector syntax must match the three-scalar syntax');

% Compare equivalent rotations through their Euler angle degeneracy
result=test_true(result,'euler_equiv zero-beta degeneracy',euler_equiv([0.3 0 0.4],[0.7 0 0],1e-14),...
                 'for beta equal to zero, alpha and gamma only enter through their sum');

% Check angular tolerance acceptance
result=test_true(result,'euler_equiv tolerance pass',euler_equiv([0 0 0],[5e-7 0 0],1e-6),...
                 'relative rotations inside the angular tolerance must pass');

% Check angular tolerance rejection
result=test_true(result,'euler_equiv tolerance fail',~euler_equiv([0 0 0],[2e-6 0 0],1e-6),...
                 'relative rotations outside the angular tolerance must fail');

% Recover a non-singular Euler rotation through its DCM
angles=[0.37 0.91 -0.42];
R=euler2dcm(angles);
angles_obs=dcm2euler(R);
result=test_close(result,'dcm2euler reconstruction',euler2dcm(angles_obs),R,1e-14,1e-14,...
                  'Euler angles are not unique, but the reconstructed DCM must be the same rotation to machine precision');
result=test_close(result,'dcm2euler beta zero',euler2dcm(dcm2euler(euler2dcm(0.4,0,1.1))),euler2dcm(0.4,0,1.1),1e-14,1e-14,...
                  'the beta=0 gimbal case must reconstruct exactly');
result=test_close(result,'dcm2euler beta pi',euler2dcm(dcm2euler(euler2dcm(0.4,pi,1.1))),euler2dcm(0.4,pi,1.1),1e-14,1e-14,...
                  'the beta=pi gimbal case must reconstruct exactly');
R_dirty=R+1e-7*[0.3 -0.7 0.2; 0.1 0.4 -0.6; -0.5 0.2 0.3];
result=test_close(result,'dcm2euler corrupted input',euler2dcm(dcm2euler(R_dirty)),R,1e-6,1e-6,...
                  'a subtly corrupted DCM must return the angles of the nearest proper rotation');

% Compose two active rotations in the documented order
rot_one=[0.2 0.4 -0.3];
rot_two=[-0.5 0.7 0.6];
rot_cmp=euler_sup(rot_one,rot_two);
R_cmp=euler2dcm(rot_cmp);
R_ref=euler2dcm(rot_two)*euler2dcm(rot_one);
result=test_close(result,'euler_sup composition',R_cmp,R_ref,1e-7,1e-7,...
                  'the composite matrix must be R_two*R_one for column vectors to optimiser tolerance');

% Check angle-axis normalisation and inverse rotation
axis_vec=[1;2;3];
angle=0.73;
R_one=anax2dcm(axis_vec,angle);
R_two=anax2dcm(5*axis_vec,angle);
R_inv=anax2dcm(axis_vec,-angle);
result=test_close(result,'anax2dcm axis normalisation',R_one,R_two,1e-15,1e-15,...
                  'multiplying the axis by a scalar must not change the represented rotation');
result=test_close(result,'anax2dcm inverse',R_one*R_inv,eye(3),1e-14,1e-14,...
                  'opposite angle-axis rotations must compose to identity');
result=test_close(result,'anax2dcm proper orthogonality',R_one'*R_one,eye(3),1e-14,1e-14,...
                  'direction cosine matrices preserve lengths');
result=test_close(result,'anax2dcm determinant',det(R_one),1,1e-14,1e-14,...
                  'proper rotations have determinant +1');
result=test_close(result,'anax2dcm active convention, Z axis',anax2dcm([0 0 1],0.31),euler2dcm(0,0,0.31),1e-14,1e-14,...
                  'a rotation about the Z axis must match the corresponding active ZYZ Euler rotation');
result=test_close(result,'anax2dcm active convention, Y axis',anax2dcm([0 1 0],0.27),euler2dcm(0,0.27,0),1e-14,1e-14,...
                  'a rotation about the Y axis must match the corresponding active ZYZ Euler rotation');

% Round-trip a non-zero angle through quaternion form
q=anax2qter(axis_vec,angle);
[axis_obs,angle_obs]=qter2anax(q);
axis_ref=(axis_vec/norm(axis_vec,2)).';
q_norm=sqrt(q.u^2+q.i^2+q.j^2+q.k^2);
result=test_close(result,'anax2qter unit norm',q_norm,1,1e-15,1e-15,...
                  'rotation quaternions must have unit norm');
result=test_close(result,'qter2anax axis',axis_obs,axis_ref,1e-15,1e-15,...
                  'for angles below pi the quaternion axis is the normalised input axis');
result=test_close(result,'qter2anax angle',angle_obs,angle,1e-15,1e-15,...
                  'for angles below pi the quaternion angle is recovered without branch ambiguity');

% Check the quaternion conversions against the Euler convention
qtr=euler2qter(0.4,1.1,-0.7);
result=test_close(result,'euler2qter unit norm',sqrt(qtr.u^2+qtr.i^2+qtr.j^2+qtr.k^2),1,1e-15,1e-15,...
                  'rotation quaternions must have unit norm');
result=test_close(result,'euler2qter/qter2dcm convention',qter2dcm(qtr),euler2dcm(0.4,1.1,-0.7),1e-14,1e-14,...
                  'the quaternion route must reproduce the active ZYZ rotation matrix');
[alp_obs,bet_obs,gam_obs]=qter2euler(qtr);
result=test_close(result,'qter2euler recovery',euler2dcm(alp_obs,bet_obs,gam_obs),euler2dcm(0.4,1.1,-0.7),1e-14,1e-14,...
                  'Euler angles are not unique, but the recovered angles must give the same rotation');
result=test_close(result,'dcm2qter/qter2dcm rotation',qter2dcm(dcm2qter(euler2dcm(0.4,1.1,-0.7))),euler2dcm(0.4,1.1,-0.7),1e-14,1e-14,...
                  'the DCM round trip through the quaternion must be exact to machine precision');

% Check that identity rotation has identity second-rank Wigner matrix
D=dcm2wigner(eye(3));
result=test_close(result,'dcm2wigner identity',D,eye(5),1e-15,1e-15,...
                  'the identity rotation leaves every second-rank tensor component unchanged');
result=test_close(result,'dcm2wigner unitarity',D'*D,eye(5),1e-15,1e-15,...
                  'Wigner rotation matrices are unitary representations of rotations');

% Align two non-zero vectors and check the proper rotation properties
v_from=[2;0;0];
v_to=[0;3;0];
R=rotmat_align(v_from,v_to);
result=test_close(result,'rotmat_align vector action',R*(v_from/norm(v_from,2)),v_to/norm(v_to,2),1e-15,1e-15,...
                  'the alignment matrix must take the normalised source vector into the normalised target vector');
result=test_close(result,'rotmat_align orthogonality',R'*R,eye(3),1e-15,1e-15,...
                  'minimum-angle alignment is a proper orthogonal rotation');
result=test_close(result,'rotmat_align determinant',det(R),1,1e-15,1e-15,...
                  'minimum-angle alignment must not introduce reflection');

% Check the anti-parallel branch explicitly
R=rotmat_align([1;0;0],[-1;0;0]);
result=test_close(result,'rotmat_align anti-parallel action',R*[1;0;0],[-1;0;0],1e-15,1e-15,...
                  'anti-parallel vectors require a pi rotation around an orthogonal axis');
result=test_close(result,'rotmat_align anti-parallel determinant',det(R),1,1e-15,1e-15,...
                  'the anti-parallel branch must still be a proper rotation');

end


