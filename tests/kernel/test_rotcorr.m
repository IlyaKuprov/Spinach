% Checks the rigid-body surface-ellipsoid diffusion estimator against
% the exact sphere limit, physical scaling, and coordinate invariants.
% Syntax:
%
%                         result=test_rotcorr()
%
% Outputs:
%
%    result - regression results and diagnostic messages
%
% talos@spindynamics.org

function result=test_rotcorr()

% Identify the physical contract
result=new_test_result('kernel/rotcorr','Surface ellipsoid diffusion',...
                      'Sphere limit, scaling, and coordinate invariants.');

% Compare with the exact rank-2 Stokes-Einstein-Debye sphere limit
[tau,D]=rotcorr([0 0 0],10,0,1.4,293,1e-3,8000);
reference=4*pi*1e-3*1e-27/(3*1.380649e-23*293);
result=test_close(result,'sphere time',tau/reference,1,0,1e-6,...
                  'the sphere has tau=4*pi*eta*r^3/(3*k*T)');
result=test_close(result,'sphere tensor',D*6*reference,eye(3),0,1e-5,...
                  'all three sphere diffusion rates coincide');

% Exercise an anisotropic overlapping-sphere body
xyz=[-3 0 0;0 0 0;3 0 0]; radii=[2;2;2];
[tau,D]=rotcorr(xyz,radii,1,1.4,300,1e-3,8000);
result=test_true(result,'positive anisotropy',...
                 all(eig(D)>0)&&(max(eig(D))/min(eig(D))>1.1),...
                 'an elongated body has positive unequal diffusion rates');

% Check rigid translation without changing the surface quadrature
[moved,tensor]=rotcorr(xyz+[17 -9 4],radii,1,1.4,300,1e-3,8000);
result=test_close(result,'translation time',moved/tau,1,0,1e-12,...
                  'translation cannot change rotational drag');
result=test_close(result,'translation tensor',tensor/norm(D),D/norm(D),...
                  0,1e-12,'translation cannot rotate the diffusion tensor');

% Check rigid rotation to within the finite surface quadrature error
angle=0.71; R=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
[moved,tensor]=rotcorr(xyz*R',radii,1,1.4,300,1e-3,8000);
result=test_close(result,'rotation time',moved/tau,1,0,3e-3,...
                  'Fibonacci quadrature converges to rotational invariance');
result=test_close(result,'rotation tensor',tensor/norm(D),R*D*R'/norm(D),...
                  0,3e-3,'the diffusion tensor transforms covariantly');

% Check temperature, viscosity, and cubic geometric scaling
scaled=rotcorr(2*xyz,2*radii,2,2.8,600,3e-3,8000);
result=test_close(result,'physical scaling',scaled/tau,12,0,1e-11,...
                  'tau scales as viscosity times length cubed over temperature');

% Check resolution convergence on a non-spherical body
fine=rotcorr(xyz,radii,1,1.4,300,1e-3,32000);
result=test_close(result,'surface convergence',fine/tau,1,0,3e-3,...
                  'refining the surface grid stabilises the scalar time');

% Reject an invalid physical solvent input
rejected=false;
try
    rotcorr(xyz,radii,1,1.4,300,0,8000);
catch
    rejected=true;
end
result=test_true(result,'zero viscosity rejected',rejected,...
                 'a positive solvent viscosity is required');

end


