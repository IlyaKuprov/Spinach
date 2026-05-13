% Tests deterministic coordinate and tensor transforms. Syntax:
%
%                    result=test_transform_roundtrip_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks transformation functions by using exact geometrical
% identities, algebraic inverses, and known tensor decompositions.
%
% ilya.kuprov@weizmann.ac.il

function result=test_transform_roundtrip_suite()

% Announce the test target
fprintf('TESTING: Coordinate and tensor transform functions\n');

% State the transform target of the test
result=new_test_result('kernel/transform_roundtrip_suite',...
                       'Coordinate and tensor transform functions',...
                       'transform helpers must preserve rotations, coordinates, and tensor decompositions.');

% Direction-cosine matrices must be orthogonal proper rotations
R=anax2dcm([0 0 2],pi/3);
result=test_close(result,'anax2dcm orthogonality',R'*R,eye(3),1e-14,1e-14,...
                  'angle-axis conversion must return an orthogonal rotation matrix');
result=test_close(result,'anax2dcm determinant',det(R),1,1e-14,1e-14,...
                  'proper rotations have determinant +1');

% Quaternion and angle-axis representations must describe the same rotation
axis=[1 2 3]; angle=0.37*pi;
q=anax2quat(axis,angle);
[axis_back,angle_back]=quat2anax(q);
R_from_axis=anax2dcm(axis,angle);
R_from_quat=anax2dcm(axis_back,angle_back);
result=test_close(result,'anax2quat/quat2anax rotation',R_from_quat,R_from_axis,1e-14,1e-14,...
                  'quaternion round-trip must preserve the represented rotation');

% Euler conversion is ill-conditioned in angles, but DCM reconstruction is unique
angles=[0.21*pi 0.37*pi 0.43*pi];
R=euler2dcm(angles);
angles_back=dcm2euler(R);
result=test_close(result,'dcm2euler/euler2dcm rotation',euler2dcm(angles_back),R,1e-7,1e-7,...
                  'Euler-angle recovery must reconstruct the original active ZYZ rotation to the documented numerical accuracy of the inverse problem');

% Axiality/rhombicity to matrix with zero Euler angles gives the Mehring-order eigenvalues
iso=4; ax=6; rh=2;
M=axrh2mat(iso,ax,rh,0,0,0);
M_ref=diag([iso-(ax+3*rh)/6, iso-(ax-3*rh)/6, iso+ax/3]);
result=test_close(result,'axrh2mat principal values',M,M_ref,1e-14,1e-14,...
                  'zero-angle axiality/rhombicity conversion must place principal values on the diagonal');
[iso_back,ax_back,rh_back,eigvals]=mat2axrh(M_ref);
result=test_close(result,'mat2axrh isotropic part',iso_back,mean(diag(M_ref)),1e-14,1e-14,...
                  'the isotropic part is the mean principal value');
result=test_close(result,'mat2axrh axiality',ax_back,eigvals(3)-(eigvals(1)+eigvals(2))/2,1e-14,1e-14,...
                  'mat2axrh axiality is zz-(xx+yy)/2 in Mehring order');
result=test_close(result,'mat2axrh rhombicity',rh_back,eigvals(1)-eigvals(2),1e-14,1e-14,...
                  'mat2axrh rhombicity is xx-yy in Mehring order');

% Cartesian and irreducible spherical tensor representations are algebraic inverses
T=[1 2 3;4 5 6;7 8 10];
[r0,r1,r2]=mat2sphten(T);
T_back=sphten2mat(r0,r1,r2);
result=test_close(result,'mat2sphten/sphten2mat',T_back,T,1e-13,1e-13,...
                  'the nine spherical tensor components span all 3x3 Cartesian tensors');

% Fractional crystallographic coordinates in an orthorhombic cell scale by cell edges
ABC=[0 0 0; 1/2 1/3 1/4; 1 1 1];
[XYZ,va,vb,vc]=frac2cart(2,3,4,90,90,90,ABC);
result=test_close(result,'frac2cart orthorhombic coordinates',XYZ,ABC*diag([2 3 4]),1e-14,1e-14,...
                  'orthorhombic fractional coordinates scale independently along a, b, and c');
result=test_close(result,'frac2cart primitive vectors',[va vb vc],diag([2 3 4]),1e-14,1e-14,...
                  'orthorhombic primitive vectors are the Cartesian unit-cell edges');

% Spherical coordinates use ISO radius/inclination/azimuth convention
[r,theta,phi]=xyz2sph([1 0 0],[0 1 0],[0 0 1]);
result=test_close(result,'xyz2sph radius',r,[1 1 1],1e-14,1e-14,...
                  'unit Cartesian basis vectors have unit radius');
result=test_close(result,'xyz2sph inclination',theta,[pi/2 pi/2 0],1e-14,1e-14,...
                  'inclination is measured down from the positive z axis');
result=test_close(result,'xyz2sph azimuth',phi,[0 pi/2 0],1e-14,1e-14,...
                  'azimuth is measured in the xy plane from the positive x axis');

end
