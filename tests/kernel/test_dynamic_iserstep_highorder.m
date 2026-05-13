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
% reference on a compact Hilbert-space problem.
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

end

% Local quiet one-spin Hilbert-space test system
function spin_system=local_hilb_system()

% Specify the spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};

% Specify a full Hilbert basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Build the quiet regression-test system
spin_system=test_spin_system(sys,inter,bas);

end


