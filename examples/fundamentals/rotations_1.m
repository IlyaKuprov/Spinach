% Tests the internal consistency of kernel rotation functions.
%
% i.kuprov@soton.ac.uk

function rotations_1()

% Generate a random symmetric traceless 3x3 matrix
A=rand(3); A=A+A'; A=A-trace(A)/3;

% Generate a random set of Euler angles
eulers=rand(1,3).*[2*pi pi 2*pi];

%% Test 1: euler2dcm, wigner, mat2sphten

% DCM rotation followed by a transformation into irreducible components
[~,~,rank2_a]=mat2sphten(euler2dcm(eulers)*A*euler2dcm(eulers)');

% Transformation into irreducible components followed by a Wigner rotation
[~,~,rank2]=mat2sphten(A); rank2_b=wigner(2,eulers(1),eulers(2),eulers(3))*rank2;

% Check the difference
if norm(rank2_a-rank2_b,2)<1e-10
    disp('Test 1 passed.');
else
    error('Test 1 inconsistency detected');
end

%% Test 2: euler2dcm, dcm2euler

% Transforms Euler angles into DCM
R=euler2dcm(eulers);

% Transform the DCM back into Euler angles
new_eulers=dcm2euler(R);

% Check the difference
if norm(eulers-new_eulers,2)<1e-3
    disp('Test 2 passed.');
else
    error('Test 2 inconsistency detected');
end

%% Test 3: euler2dcm, wigner, dcm2wigner

% Transforms Euler angles into DCM, then DCM to Wigner matrix
W_a=dcm2wigner(euler2dcm(eulers));

% Transform Euler angles to Wigner matrix
W_b=wigner(2,eulers(1),eulers(2),eulers(3));

% Check the difference
if norm(W_a-W_b,2)<1e-10
    disp('Test 3 passed.');
else
    error('Test 3 inconsistency detected');
end

%% Test 4: mat2sphten, sphten2mat

% Transform the matrix into spherical tensors
[rank0,rank1,rank2]=mat2sphten(A);

% Transform back
B=sphten2mat(rank0,rank1,rank2);

% Check the difference
if norm(A-B,2)<1e-10
    disp('Test 4 passed.');
else
    error('Test 4 inconsistency detected');
end

%% Test 5: dcm2quat, quat2dcm

% Generate a DCM
R=euler2dcm(eulers);

% Convert to quaternions and back
if norm(R-quat2dcm(dcm2quat(R)),2)<1e-10
    disp('Test 5 passed.');
else
    error('Test 5 inconsistency detected');
end

%% Test 6: anax2quat, quat2anax

% Generate a unit quaternion
q1.u=randn(); q1.i=randn();
q1.j=randn(); q1.k=randn();
qnorm=norm([q1.u q1.i q1.j q1.k],2);
q1.u=q1.u/qnorm; q1.i=q1.i/qnorm;
q1.j=q1.j/qnorm; q1.k=q1.k/qnorm;

% Convert to angle-axis and back
[aa_axis,aa_angle]=quat2anax(q1);
q2=anax2quat(aa_axis,aa_angle);

% Check the difference
if norm([q1.u-q2.u q1.i-q2.i q1.j-q2.j q1.k-q2.k],2)<1e-10
    disp('Test 6 passed.');
else
    error('Test 6 inconsistency detected');
end

%% Test 7: quat2anax, anax2dcm, quat2dcm

% Generate a unit quaternion
q.u=randn(); q.i=randn();
q.j=randn(); q.k=randn();
qnorm=norm([q.u q.i q.j q.k],2);
q.u=q.u/qnorm; q.i=q.i/qnorm;
q.j=q.j/qnorm; q.k=q.k/qnorm;

% Convert directly to DCM
R1=quat2dcm([q.u q.i q.j q.k]);

% Convert to DCM via angle-axis
[aa_axis,aa_angle]=quat2anax(q);
R2=anax2dcm(aa_axis,aa_angle);

% Check the difference
if norm(R1-R2,2)<1e-10
    disp('Test 7 passed.');
else
    error('Test 7 inconsistency detected');
end

end

