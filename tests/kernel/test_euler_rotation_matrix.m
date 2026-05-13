% Tests active ZYZ Euler rotation matrices. Syntax:
%
%                    result=test_euler_rotation_matrix()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the active convention used by Spinach: alpha=pi/2,
% beta=0, gamma=0 is a counter-clockwise rotation around Z, taking x into y.
%
% ilya.kuprov@weizmann.ac.il

function result=test_euler_rotation_matrix()

% Announce the test target
fprintf('TESTING: Euler active rotation matrix\n');

% State the physical target of the test
result=new_test_result('kernel/euler_rotation_matrix',...
                       'Euler active rotation matrix',...
                       'euler2dcm must implement the active ZYZ convention.');

% Build a simple ninety-degree Z rotation
R=euler2dcm(pi/2,0,0);
R_ref=[0 -1 0;1 0 0;0 0 1];

% Check the active rotation and orthogonality
result=test_close(result,'active Z rotation',R,R_ref,1e-15,1e-15,...
                  'a positive active Z rotation maps the x axis into y');
result=test_close(result,'orthogonality',R'*R,eye(3),1e-15,1e-15,...
                  'proper rotations preserve vector lengths');
result=test_close(result,'proper determinant',det(R),1,1e-15,1e-15,...
                  'direction cosine matrices must have determinant +1');
result=test_close(result,'vector action',R*[1;0;0],[0;1;0],1e-15,1e-15,...
                  'the documented action is v=R*v for column vectors');

end

