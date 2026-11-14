% Tests nonlinear high-order iserstep branches. Syntax:
%
%               result=test_dynamic_iserstep_highorder()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks zero-step handling, nonlinear generator execution, and
% agreement of high-order Lie and RKMK branches against a refined DP8
% reference on a compact Hilbert-space problem. Additional LG4A checks
% cover constant and commuting fields, fourth-order convergence in a
% rotating field, and a genuinely nonlinear density-matrix flow against
% independently refined ODE references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_iserstep_highorder()

% Announce the test target
fprintf('TESTING: Nonlinear high-order iserstep branches\n');

% State the nonlinear Lie-step target of the test
result=new_test_result('kernel/dynamic_iserstep_highorder',...
                       'Nonlinear high-order iserstep branches',...
                       'high-order Lie and RKMK solvers must execute nonlinear generators while preserving density-matrix invariants.');

% Build a one-proton Hilbert-space spin system
spin_system=local_hilb_system();
S=pauli(2);
Sx=S.x;
Sy=S.y;
Sz=S.z;

% Build a Hermitian density matrix with non-zero coherences
rho=[0.7,0.2+0.1i;0.2-0.1i,0.3];

% Define a mildly nonlinear, non-commuting Hamiltonian field
H0=0.31*Sx+0.17*Sz;
H1=-0.23*Sy+0.05*Sx;
H2=0.11*Sx-0.07*Sz;
Lfun=@(t,state)H0+sin(4*t)*H1+0.2*real(trace(Sz*state))*H2;
dt=5e-3;

% Check the explicit zero-time shortcut in LG4A
rho_zero=iserstep(spin_system,{Lfun,0,'LG4A'},rho,0);
result=test_close(result,'iserstep LG4A zero step',rho_zero,rho,1e-15,1e-15,...
                  'the LG4A branch must return the input state for a zero time step');

% Build a refined DP8 reference using two half steps
rho_ref=iserstep(spin_system,{Lfun,0,'RKMK-DP8'},rho,dt/2);
rho_ref=iserstep(spin_system,{Lfun,dt/2,'RKMK-DP8'},rho_ref,dt/2);

% Exercise nonlinear high-order branches against the refined reference
methods={'LG4A','RKMK4','RKMK-DP5','RKMK-DP8','RKMK-RKF45'};
tolerances=[5e-8 5e-8 5e-9 5e-10 5e-9];
for n=1:numel(methods)

    % Take one nonlinear Lie step with this method
    rho_obs=iserstep(spin_system,{Lfun,0,methods{n}},rho,dt);

    % Compare to the refined high-order reference
    result=test_close(result,['iserstep nonlinear ' methods{n}],rho_obs,rho_ref,tolerances(n),tolerances(n),...
                      'nonlinear high-order branches agree with a refined DP8 reference on a short time step');

    % Check density-matrix invariants preserved by unitary propagation
    result=test_close(result,['iserstep trace ' methods{n}],trace(rho_obs),trace(rho),1e-12,1e-12,...
                      'unitary Hamiltonian propagation preserves the density-matrix trace');
    result=test_close(result,['iserstep Hermiticity ' methods{n}],rho_obs,rho_obs',1e-12,1e-12,...
                      'unitary Hamiltonian propagation preserves Hermiticity');
end

% Check the constant-generator limit against exact density conjugation
constant_fun=@(t,state)H0;
rho_const=iserstep(spin_system,{constant_fun,0,'LG4A'},rho,0.2);
P=expm(-1i*0.2*H0);
result=test_close(result,'LG4A constant generator',rho_const,P*rho*P',1e-13,1e-13,...
                  'a constant Hamiltonian must give exact density conjugation');

% Check exact Simpson phase and trace for a commuting quadratic field
quadratic_fun=@(t,state)(1+t^2)*Sz;
rho_obs=rho;
for n=1:8
    rho_obs=iserstep(spin_system,{quadratic_fun,(n-1)/8,'LG4A'},rho_obs,1/8);
end
P=expm(-1i*(4/3)*Sz);
result=test_close(result,'LG4A commuting quadratic',rho_obs,P*rho*P',1e-12,1e-12,...
                  'the commuting quadratic field has an exact Simpson phase and no artificial damping');

% Specify a rotating noncommuting field with an exact frame solution
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=diag([1 -1]);
rotating_fun=@(t,state)0.8*(cos(2.3*t)*X+sin(2.3*t)*Y)+0.4*Z;
P=expm(-1i*2.3*Z)*expm(-1i*2*(0.8*X+(0.4-2.3/2)*Z));
rho_exact=P*rho*P';
grid_sizes=[16 32 64];
rotating_errors=zeros(size(grid_sizes));

% Refine the same physical interval without projecting the density matrix
for n=1:numel(grid_sizes)

    % Integrate the rotating field over two seconds
    time_step=2/grid_sizes(n);
    rho_obs=rho;
    for k=1:grid_sizes(n)
        rho_obs=iserstep(spin_system,{rotating_fun,(k-1)*time_step,'LG4A'},rho_obs,time_step);
    end
    rotating_errors(n)=norm(rho_obs-rho_exact,'fro');

end
rotating_orders=log2(rotating_errors(1:end-1)./rotating_errors(2:end));
result=test_close(result,'LG4A rotating field order',rotating_orders,4*ones(1,2),0.05,0,...
                  'a noncommuting prescribed field must exhibit fourth-order global convergence');

% Build independent references for a genuinely state-dependent Hamiltonian
nonlinear_fun=@(t,state)0.7*X+0.2*cos(0.9*t)*Y+...
                        (0.3+0.6*real(trace(Z*state)))*Z;
loose_opts=odeset('RelTol',1e-11,'AbsTol',1e-12);
tight_opts=odeset('RelTol',2.5e-14,'AbsTol',2.5e-15);
[~,ref_loose]=ode113(@nonlin_density_rhs,[0 1.5],rho(:),loose_opts);
[~,ref_tight]=ode113(@nonlin_density_rhs,[0 1.5],rho(:),tight_opts);
[~,ref_other]=ode89(@nonlin_density_rhs,[0 1.5],rho(:),tight_opts);
rho_ref=reshape(ref_tight(end,:).',2,2);
result=test_close(result,'LG4A nonlinear reference refinement',ref_loose(end,:),ref_tight(end,:),2e-10,0,...
                  'tightening the independent ODE reference must leave the final density matrix unchanged');
result=test_close(result,'LG4A independent nonlinear reference',ref_other(end,:),ref_tight(end,:),2e-12,0,...
                  'independent ODE methods must agree well below the measured LG4A errors');
nonlinear_errors=zeros(size(grid_sizes));

% Check full-interval nonlinear order and physical density invariants
for n=1:numel(grid_sizes)

    % Evaluate state feedback at every native LG4A stage
    time_step=1.5/grid_sizes(n);
    rho_obs=rho;
    for k=1:grid_sizes(n)
        rho_obs=iserstep(spin_system,{nonlinear_fun,(k-1)*time_step,'LG4A'},rho_obs,time_step);
    end
    nonlinear_errors(n)=norm(rho_obs-rho_ref,'fro');
    result=test_close(result,['LG4A nonlinear trace ' num2str(grid_sizes(n))],...
                      trace(rho_obs),trace(rho),1e-12,1e-12,...
                      'the nonlinear flow must preserve density-matrix trace');
    result=test_close(result,['LG4A nonlinear spectrum ' num2str(grid_sizes(n))],...
                      sort(eig(rho_obs)),sort(eig(rho)),1e-12,1e-12,...
                      'unitary conjugation must preserve the density-matrix eigenvalues');

end
nonlinear_orders=log2(nonlinear_errors(1:end-1)./nonlinear_errors(2:end));
result=test_close(result,'LG4A nonlinear field order',nonlinear_orders,4*ones(1,2),0.08,0,...
                  'a genuinely state-dependent Hamiltonian must exhibit fourth-order global convergence');

end

% Independent density-matrix reference equation
function derivative=nonlin_density_rhs(time,rho_vector)

% Reconstruct the Hermitian mean-field Hamiltonian
rho=reshape(rho_vector,2,2);
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=diag([1 -1]);
H=0.7*X+0.2*cos(0.9*time)*Y+(0.3+0.6*real(trace(Z*rho)))*Z;

% Integrate the commutator directly without exponentiating
derivative=reshape(-1i*(H*rho-rho*H),4,1);

end

% Local quiet one-spin Hilbert-space test system
function spin_system=local_hilb_system()

% Specify the spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};

% Specify a full Hilbert basis
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Build the quiet regression-test system
spin_system=test_spin_system(sys,inter,bas);

end


