% Tests tensor transform helpers. Syntax:
%
%                    result=test_transform_tensor_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks interaction tensor parametrisations, spherical tensor
% round-trips, quadrupolar conversions, axial symmetrisation, and simple
% Hamiltonian decomposition.
%
% ilya.kuprov@weizmann.ac.il

function result=test_transform_tensor_suite()

% Announce the test target
fprintf('TESTING: Tensor transform helpers\n');

% State the tensor target of the test
result=new_test_result('kernel/transform_tensor_suite',...
                       'Tensor transform helpers',...
                       'tensor transforms must preserve their algebraic definitions and round-trips.');

% Check Haeberlen anisotropy and asymmetry at zero Euler angles
iso=10; aniso=6; asym=0.25;
red_aniso=2*aniso/3;
M_ref=diag([iso-red_aniso*(1+asym)/2,iso-red_aniso*(1-asym)/2,iso+red_aniso]);
M=anas2mat(iso,aniso,asym,0,0,0);
result=test_close(result,'anas2mat principal values',M,M_ref,1e-15,1e-15,...
                  'zero Euler angles must return the Haeberlen principal-axis tensor');

% Check axiality and rhombicity matrix construction
iso=4; ax=6; rh=2;
M_ref=diag([2 4 6]);
M=axrh2mat(iso,ax,rh,0,0,0);
result=test_close(result,'axrh2mat principal values',M,M_ref,1e-15,1e-15,...
                  'axiality 2*zz-(xx+yy) and rhombicity yy-xx define the principal values');
[iso_obs,ax_obs,rh_obs,eigs_obs]=mat2axrh(M_ref);
result=test_close(result,'mat2axrh isotropic part',iso_obs,4,1e-15,1e-15,...
                  'the isotropic part is the average of the three eigenvalues');
result=test_close(result,'mat2axrh axiality',ax_obs,6,1e-15,1e-15,...
                  'mat2axrh reports 2*zz-(xx+yy) for Mehring-ordered eigenvalues');
result=test_close(result,'mat2axrh rhombicity',rh_obs,2,1e-15,1e-15,...
                  'mat2axrh reports yy-xx for Mehring-ordered eigenvalues');
result=test_close(result,'mat2axrh eigenvalue order',eigs_obs,[2;4;6],1e-15,1e-15,...
                  'Mehring order sorts the principal values as xx<=yy<=zz');

% Check span and skew construction in the Herzfeld-Berger convention
iso=5; span=6; skew=0.5;
M_ref=diag([1.5 6 7.5]);
M=spsk2mat(iso,span,skew,0,0,0);
result=test_close(result,'spsk2mat principal values',M,M_ref,1e-15,1e-15,...
                  'span is zz-xx, and skew fixes the middle eigenvalue displacement from isotropy');

% Check zero-field splitting tensor construction
D=9; E=2;
M_ref=diag([-1 -5 6]);
M=zfs2mat(D,E,0,0,0);
result=test_close(result,'zfs2mat principal values',M,M_ref,1e-15,1e-15,...
                  'D and E define [-D/3+E,-D/3-E,2D/3] in the principal frame');
result=test_close(result,'zfs2mat tracelessness',trace(M),0,1e-15,1e-15,...
                  'zero-field splitting tensors are traceless by definition');

% Check isotropic-antisymmetric-symmetric decomposition and reconstruction
C=[1 2 -3;4 -5 6;7 8 9];
[a,d,A]=mat2ias(C);
C_obs=ias2mat(a,d,A);
result=test_close(result,'mat2ias isotropic scalar',a,trace(C)/3,1e-15,1e-15,...
                  'the isotropic scalar is one third of the Cartesian trace');
result=test_close(result,'mat2ias symmetric traceless part',trace(A),0,1e-15,1e-15,...
                  'the symmetric residual must be traceless after removing the isotropic part');
result=test_close(result,'mat2ias symmetric part',A,A.',1e-15,1e-15,...
                  'the rank-two residual of the IAS decomposition is symmetric');
result=test_close(result,'ias2mat round-trip',C_obs,C,1e-15,1e-15,...
                  'IAS reconstruction must recover every element of the original 3x3 matrix');

% Round-trip the Cartesian tensor through irreducible spherical components
[rank0,rank1,rank2]=mat2sphten(C);
C_obs=sphten2mat(rank0,rank1,rank2);
result=test_close(result,'mat2sphten sphten2mat round-trip',C_obs,C,1e-15,1e-15,...
                  'the nine spherical tensor coefficients form a complete Cartesian tensor basis');

% Check spherical harmonic coefficients of an isotropic quadratic form
[r0,r1,r2]=qform2sph(3*eye(3));
result=test_close(result,'qform2sph isotropic rank zero',r0,6*sqrt(pi),1e-15,1e-15,...
                  'an isotropic quadratic form contributes only to Y00 with coefficient 2*sqrt(pi)*a');
result=test_close(result,'qform2sph isotropic rank one',r1,zeros(1,3),1e-15,1e-15,...
                  'a symmetric quadratic form has no vector spherical harmonic component');
result=test_close(result,'qform2sph isotropic rank two',r2,zeros(1,5),1e-15,1e-15,...
                  'an isotropic quadratic form has no rank-two anisotropy');

% Check the axial Stevens q=0 component through the rank-one transform
Bkq=stev2sph(1,[0;1;0]);
result=test_close(result,'stev2sph rank-one axial component',Bkq,[0;1;0],1e-15,1e-15,...
                  'the q=0 rank-one Stevens component is already the axial spherical component');

% Check traceless symmetric matrix parameter extraction
T=diag([-2 -1 3]);
[ax_obs,rh_obs,angles_obs]=tsm2param(T);
result=test_close(result,'tsm2param axiality',ax_obs,9,1e-12,1e-12,...
                  'Mehring axiality for [-2,-1,3] is 2*3-(-2-1)');
result=test_close(result,'tsm2param rhombicity',rh_obs,1,1e-12,1e-12,...
                  'Mehring rhombicity for [-2,-1,3] is yy-xx');
result=test_close(result,'tsm2param orientation reconstruction',euler2dcm(angles_obs)*T*euler2dcm(angles_obs).',T,1e-7,1e-7,...
                  'the returned Euler angles must reconstruct the principal-frame tensor to optimiser tolerance');

% Check quadrupolar tensor construction at zero Euler angles
Cq=1200; eta=0.2; spin_q=1;
Q_ref=diag([-240 -360 600]);
Q=eeqq2nqi(Cq,eta,spin_q,[0 0 0]);
result=test_close(result,'eeqq2nqi principal values',Q,Q_ref,1e-12,1e-15,...
                  'Cq and eta define a traceless quadrupolar tensor with principal values XX, YY, and ZZ');

% Check CASTEP electric-field-gradient scaling convention
V=diag([-1 -2 3]);
quad_moment=0.2;
scale=9.717362e+21*(quad_moment*1e-28)*1.60217657e-19/(6.62606957e-34*2*spin_q*(2*spin_q-1));
Q=castep2nqi(V,quad_moment,spin_q);
result=test_close(result,'castep2nqi scaling',Q,scale*V,1e-9,1e-15,...
                  'CASTEP atomic-unit EFG tensors are scaled by e*q*Q/h/[2I(2I-1)]');

% Check WebLab cone convention delegation for two sites
alpha=0.11; theta=0.42; phi=0.73;
[Q1,Q2]=weblab2nqi(Cq,eta,spin_q,alpha,theta,phi);
result=test_close(result,'weblab2nqi first two-site tensor',Q1,eeqq2nqi(Cq,eta,spin_q,[-phi/2 theta alpha]),1e-12,1e-15,...
                  'the first two-site WebLab tensor uses Euler angles [-phi/2 theta alpha]');
result=test_close(result,'weblab2nqi second two-site tensor',Q2,eeqq2nqi(Cq,eta,spin_q,[+phi/2 theta alpha]),1e-12,1e-15,...
                  'the second two-site WebLab tensor uses Euler angles [+phi/2 theta alpha]');

% Check rotational averaging around the z axis
T=diag([1 3 5]);
T_obs=axis_tsymm(T,[0;0;1]);
T_ref=diag([2 2 5]);
result=test_close(result,'axis_tsymm z-axis average',T_obs,T_ref,1e-12,1e-12,...
                  'full rotation around z averages xx and yy while preserving zz');

% Check spin-half Hamiltonian decomposition into Zeeman terms
S=pauli(2);
omega=[11 -7 5];
H=full(omega(1)*S.x+omega(2)*S.y+omega(3)*S.z);
[omega_obs,Q_obs]=ham2nqi(H);
result=test_close(result,'ham2nqi spin-half Zeeman vector',omega_obs,omega,1e-12,1e-12,...
                  'a spin-half traceless Hermitian Hamiltonian is fully described by its Zeeman vector');
result=test_close(result,'ham2nqi spin-half quadrupole zero',Q_obs,zeros(3),1e-15,1e-15,...
                  'spin one half has no independent quadrupolar tensor term');

end

