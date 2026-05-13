% Tests dynamic relaxation model helper paths. Syntax:
%
%              result=test_dynamic_relaxation_models_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks anisotropic and functional T1/T2 rates, damping, Lindblad,
% scalar Redfield, correlation functions, and the serial Redfield include.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_relaxation_models_suite()

% Announce the test target
fprintf('TESTING: Dynamic relaxation model helpers\n');

% State the relaxation-model target of the test
result=new_test_result('kernel/dynamic_relaxation_models_suite',...
                       'Dynamic relaxation model helpers',...
                       'relaxation helpers must assign mathematically expected rates on tiny spin systems.');

% Check anisotropic tensor rates in the extended T1/T2 model
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
inter.relaxation={'t1_t2'};
inter.r1_rates={diag([1 3 5])};
inter.r2_rates={diag([2 4 6])};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.temperature=300;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_t=test_spin_system(sys,inter,bas);
[r1_op,r2_op]=rlx_t1_t2(spin_t,[0 0 0]);
rho_z=state(spin_t,'Lz','1H');
rho_p=state(spin_t,'L+','1H');
result=test_close(result,'anisotropic T1 tensor branch',r1_op*rho_z,-5*rho_z,1e-12,1e-12,...
                  'at zero Euler angles the laboratory z axis samples the zz element of the R1 tensor');
result=test_close(result,'anisotropic T2 tensor branch',r2_op*rho_p,-6*rho_p,1e-12,1e-12,...
                  'at zero Euler angles the laboratory z axis samples the zz element of the R2 tensor');

% Check function-handle rates in the extended T1/T2 model
inter.r1_rates={@(alp,bet,gam)2+0*alp+0*bet+0*gam};
inter.r2_rates={@(alp,bet,gam)3+0*alp+0*bet+0*gam};
spin_f=test_spin_system(sys,inter,bas);
[r1_op,r2_op]=rlx_t1_t2(spin_f,[0.2 0.3 0.4]);
rho_z=state(spin_f,'Lz','1H');
rho_p=state(spin_f,'L+','1H');
result=test_close(result,'functional T1 branch',r1_op*rho_z,-2*rho_z,1e-12,1e-12,...
                  'a constant periodic function handle must set the R1 rate directly');
result=test_close(result,'functional T2 branch',r2_op*rho_p,-3*rho_p,1e-12,1e-12,...
                  'a constant periodic function handle must set the R2 rate directly');

% Check non-selective damping in Liouville space
inter_d.zeeman.scalar={0};
inter_d.relaxation={'damp'};
inter_d.damp_rate=5;
inter_d.equilibrium='zero';
inter_d.rlx_keep='labframe';
inter_d.temperature=300;
spin_d=test_spin_system(sys,inter_d,bas);
R=relaxation(spin_d);
rho_u=unit_state(spin_d);
rho_p=state(spin_d,'L+','1H');
result=test_close(result,'damping leaves unit state',R*rho_u,0*rho_u,1e-12,1e-12,...
                  'Liouville-space damping must not damp the thermodynamic unit state');
result=test_close(result,'damping transverse rate',R*rho_p,-5*rho_p,1e-12,1e-12,...
                  'non-selective damping assigns the same negative rate to non-unit states');

% Check one-spin Lindblad relaxation rates
inter_l.zeeman.scalar={0};
inter_l.relaxation={'lindblad'};
inter_l.lind_r1_rates=4;
inter_l.lind_r2_rates=7;
inter_l.equilibrium='zero';
inter_l.rlx_keep='labframe';
inter_l.temperature=300;
spin_l=test_spin_system(sys,inter_l,bas);
R=relaxation(spin_l);
rho_u=unit_state(spin_l);
rho_z=state(spin_l,'Lz','1H');
rho_p=state(spin_l,'L+','1H');
result=test_close(result,'Lindblad leaves unit state',R*rho_u,0*rho_u,1e-12,1e-12,...
                  'a trace-preserving Lindblad generator must leave the unit state unchanged');
result=test_close(result,'Lindblad longitudinal rate',R*rho_z,-4*rho_z,1e-12,1e-12,...
                  'the one-spin Lindblad R1 rate damps longitudinal magnetisation');
result=test_close(result,'Lindblad transverse rate',R*rho_p,-7*rho_p,1e-12,1e-12,...
                  'the one-spin Lindblad R2 rate damps transverse magnetisation');

% Check scalar Redfield integral against the closed zero-H0 reference
spin_l.tols.rlx_integration=1e-5;
H0=sparse(2,2);
H1=sparse([0 1; 1 0]);
tau_c=1e-3;
R=rlx_scalar(spin_l,H0,H1,{[2 tau_c]});
scalar_ref=-2*tau_c*(1-spin_l.tols.rlx_integration^2)*speye(2);
result=test_close(result,'scalar Redfield zero-H0 integral',R,scalar_ref,1e-11,1e-11,...
                  'with H0=0 and H1^2=1 the scalar Redfield integral is analytic');

% Build a Redfield-ready spin system for correlation-function checks
sys_r.magnet=14.1;
sys_r.isotopes={'1H'};
sys_r.disable={'hygiene','asyredf'};
inter_r.zeeman.scalar={0};
inter_r.relaxation={'redfield'};
inter_r.tau_c={1e-9};
inter_r.equilibrium='zero';
inter_r.rlx_keep='labframe';
inter_r.rlx_dfs='ignore';
inter_r.temperature=300;
spin_r=test_spin_system(sys_r,inter_r,bas);

% Check isotropic rank-two correlation-function weights and rates
[weights,rates,states]=corrfun(spin_r,2,3,2,3,2);
result=test_close(result,'isotropic corrfun weight',weights{1},1/5,1e-15,1e-15,...
                  'isotropic rank-two Wigner autocorrelation has weight 1/(2L+1)');
result=test_close(result,'isotropic corrfun rate',rates{1},-1e9,1e-6,1e-15,...
                  'a supplied second-rank isotropic correlation time gives rate -1/tau_c');
result=test_close(result,'corrfun species state count',nnz(states{1}),3,0,0,...
                  'a one-spin spherical-tensor basis has three non-unit states in its chemical species');

% Check axial rotational diffusion rate formula
spin_r.rlx.tau_c={[2e-9 5e-9]};
[weights,rates]=corrfun(spin_r,2,3,2,3,2);
D_ax=1/(6*2e-9);
D_eq=1/(6*5e-9);
rate_ref=-(6*D_eq+(2-2+1)^2*(D_ax-D_eq));
result=test_close(result,'axial corrfun weight',weights{1},1/5,1e-15,1e-15,...
                  'axial diffusion keeps the diagonal Wigner autocorrelation weight 1/(2L+1)');
result=test_close(result,'axial corrfun rate',rates{1},rate_ref,1e-6,1e-15,...
                  'axial diffusion rates follow the analytical rank-two expression');

% Check rhombic rotational diffusion rate formula
spin_r.rlx.tau_c={[1e-9 2e-9 3e-9]};
[~,rates]=corrfun(spin_r,2,3,3,3,3);
Dxx=1/(6e-9);
Dyy=1/(12e-9);
Dzz=1/(18e-9);
delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz);
rates_ref=-( [4*Dxx+Dyy+Dzz, Dxx+4*Dyy+Dzz, Dxx+Dyy+4*Dzz, ...
              2*Dxx+2*Dyy+2*Dzz-2*delta, 2*Dxx+2*Dyy+2*Dzz+2*delta] );
result=test_close(result,'rhombic corrfun rates',rates{1},rates_ref,1e-6,1e-15,...
                  'rhombic rank-two rotational diffusion returns the five analytical rates');

% Exercise the serial Redfield include on a tiny anisotropic one-spin system
inter_r.zeeman.matrix={diag([1 2 -3])};
spin_red=test_spin_system(sys_r,inter_r,bas);
R=relaxation(spin_red);
rho_u=unit_state(spin_red);
result=test_true(result,'serial Redfield nonzero matrix',nnz(R)>0,...
                 'anisotropic one-spin Redfield relaxation must produce nonzero matrix elements');
result=test_close(result,'serial Redfield trace preservation',R*rho_u,0*rho_u,1e-12,1e-12,...
                  'the Redfield relaxation generator must preserve the thermodynamic unit state');

end


