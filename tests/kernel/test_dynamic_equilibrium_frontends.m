% Tests equilibrium and residual-order dynamic front-end kernels. Syntax:
%
%                    result=test_dynamic_equilibrium_frontends()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises thermalize(), steady(), and residual() on compact
% Liouville-space systems with explicit fixed-point references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_equilibrium_frontends()

% State the dynamic equilibrium target of the test
result=new_test_result('kernel/dynamic_equilibrium_frontends',...
                       'Dynamic equilibrium front ends',...
                       'thermalisation and steady-state helpers must produce explicit fixed points.');

% Check IME and DiBari thermalisation branches
result=local_test_thermalize(result);

% Check Newton and squaring steady-state solvers
result=local_test_steady(result);

% Check weak residual-order tensor reduction
result=local_test_residual(result);

end


function result=local_test_thermalize(result)

% Build a one-spin spherical-tensor Liouville-space system
spin_system=local_liouville_system();
dim=size(spin_system.bas.basis,1);
unit=unit_state(spin_system);
R=-diag([0;ones(dim-1,1)]);
R=sparse(R);
rho_eq=unit+0.1*state(spin_system,'Lz','1H');

% Thermalise by the inhomogeneous master equation route
R_ime=thermalize(spin_system,R,[],[],rho_eq,'IME');
result=test_close(result,'thermalize IME fixed point',R_ime*rho_eq,zeros(dim,1),1e-14,1e-14,...
                  'thermalize() IME mode must make the requested equilibrium state stationary');

% Thermalise by the DiBari-Levitt route and compare to the defining product
H_left=operator(spin_system,'Lz','1H','left');
temperature=300.0;
R_dibari=thermalize(spin_system,R,H_left,temperature,[],'dibari');
beta=spin_system.tols.hbar/(spin_system.tols.kbol*temperature);
R_ref=R*propagator(spin_system,H_left,1i*beta);
result=test_close(result,'thermalize DiBari product',R_dibari,R_ref,1e-14,1e-14,...
                  'thermalize() DiBari mode must multiply relaxation by the imaginary-time propagator');

end


function result=local_test_steady(result)

% Build a one-spin spherical-tensor Liouville-space system
spin_system=local_liouville_system();
dim=size(spin_system.bas.basis,1);

% Construct a contractive affine propagator with a known fixed point
rho_ss=zeros(dim,1);
rho_ss(1)=1.0;
rho_ss(2)=0.20;
rho_ss(3)=-0.10;
rho_ss(4)=0.05;
contract=diag([0.25 0.50 0.75]);
P=zeros(dim);
P(1,1)=1.0;
P(2:end,1)=(eye(dim-1)-contract)*rho_ss(2:end);
P(2:end,2:end)=contract;
P=sparse(P);

% Check that the constructed propagator has the intended steady state
result=test_close(result,'steady constructed fixed point',P*rho_ss,rho_ss,1e-14,1e-14,...
                  'the test propagator is contractive around the reference fixed point');

% Recover the fixed point by Newton iteration
rho_obs=steady(spin_system,P,[],'newton');
result=test_close(result,'steady Newton solver',rho_obs,rho_ss,1e-12,1e-12,...
                  'steady() Newton mode must recover the fixed point of the propagator');

% Recover the fixed point by repeated squaring
rho_obs=steady(spin_system,P,[],'squaring');
result=test_close(result,'steady squaring solver',rho_obs,rho_ss,1e-10,1e-10,...
                  'steady() squaring mode must converge to the same fixed point');

end


function result=local_test_residual(result)

% Build a heteronuclear spin system with coordinates and weak order
sys.magnet=5.9;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={5.0,65.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=140.0;
inter.coupling.scalar{2,2}=0.0;
inter.coordinates={[0.0 0.0 0.0];[0.6 0.7 0.8]};
inter.order_matrix={diag([1e-3 2e-3 -3e-3])};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Store the original coupling tensor and apply residual ordering
J_before=spin_system.inter.coupling.matrix{1,2};
spin_system=residual(spin_system);
J_after=spin_system.inter.coupling.matrix{1,2};

% Check isotropic trace preservation and axial residual tensor shape
result=test_close(result,'residual coupling trace',trace(J_after),trace(J_before),1e-10,1e-12,...
                  'residual() must preserve the isotropic part of each coupling tensor');
result=test_close(result,'residual coupling xy degeneracy',J_after(1,1),J_after(2,2),1e-12,1e-12,...
                  'residual() must produce an axially symmetric residual tensor');
result=test_close(result,'residual coupling off-diagonal removal',J_after-diag(diag(J_after)),zeros(3),1e-12,1e-12,...
                  'residual() must remove off-diagonal tensor components after weak ordering');

end


function spin_system=local_liouville_system()

% Build the common one-spin spherical-tensor Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0.0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'nmr');

end


