% Tests the Markovian back-projection of stochastic Liouville equation
% generators onto the spin Liouville space. Syntax:
%
%                        result=test_sle2rlx()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks sle2rlx() against two independently known answers:
% the fast-exchange linewidth of two-site chemical exchange, with the
% exact eigenvalue of the exchange generator as the reference, and the
% Redfield relaxation superoperator of a single spin with an anisotro-
% pic chemical shift, with the gridfree() SLE generator as the input.
%
% talos@spindynamics.org
function result=test_sle2rlx()

% Announce the test target
fprintf('TESTING: SLE back-projection onto spin space\n');

% Start a result record
result=new_test_result('kernel/sle2rlx','SLE back-projection onto spin space',...
                       'sle2rlx must reproduce exact exchange lines and Redfield theory');

% Two-site exchange of a single coherence
omega_a=2*pi*100; omega_b=2*pi*300; k_ab=4000; k_ba=6000;
L=sparse(diag([omega_a omega_b])+1i*[-k_ab k_ba; k_ab -k_ba]);
p_a=k_ba/(k_ab+k_ba); p_b=k_ab/(k_ab+k_ba);

% Minimal spin system for the reporting interface
sys.isotopes={'1H'}; sys.magnet=1; sys.output='hush';
inter=struct(); bas.formalism='sphten-liouv'; bas.approximation={'none'};
spin_system=create(sys,inter); spin_system=basis(spin_system,bas);

% Back-projection onto the species-averaged coherence
[R,L0]=sle2rlx(spin_system,L,[p_a; p_b],[1; 1]);

% Check the species-averaged frequency
result=test_close(result,'exchange averaged frequency',L0,...
                  p_a*omega_a+p_b*omega_b,0,1e-12,...
                  'the lattice-averaged Liouvillian must be the population-weighted frequency');

% Check the fast-exchange linewidth against the exact eigenvalue
spin_line=eig(full(L)); [~,idx]=min(abs(spin_line-L0)); spin_line=spin_line(idx);
delta_omega=omega_b-omega_a; k_ex=k_ab+k_ba;
result=test_close(result,'exchange linewidth',L0+1i*R,spin_line,...
                  0,(delta_omega/k_ex)^2,...
                  'the Markov generator must match the exact spin line to the second order in the exchange parameter');

% Check the Swift-Connick limit of the exchange broadening
result=test_close(result,'exchange Swift-Connick rate',real(R),...
                  -p_a*p_b*delta_omega^2/k_ex,0,(delta_omega/k_ex)^2,...
                  'the fast-exchange broadening must match the Swift-Connick formula');

% Single spin with an anisotropic chemical shift
clear('sys','inter');
sys.magnet=14.1; sys.isotopes={'13C'}; sys.output='hush';
inter.zeeman.matrix={[200 0 0; 0 50 0; 0 0 -70]};
inter.relaxation={'redfield'}; inter.equilibrium='zero';
inter.rlx_keep='labframe'; inter.rlx_dfs='keep'; inter.tau_c={1e-11};
spin_system=create(sys,inter); spin_system=basis(spin_system,bas);

% Redfield superoperator
R_rf=relaxation(spin_system);

% SLE generator in the lab frame
parameters.spins={'13C'}; parameters.max_rank=4; parameters.tau_c=1e-11;
parameters.verbose=0;
capture=@(spin_system,parameters,H,R,K)H+1i*R+1i*K;
L=gridfree(spin_system,capture,parameters,'labframe');

% Back-projection onto the rank zero Wigner function
P=sparse(1,1,1,size(L,1)/size(R_rf,1),1);
[R_sle,L0]=sle2rlx(spin_system,L,P,P);

% Check the lattice-averaged Liouvillian
H0=hamiltonian(assume(spin_system,'labframe'));
result=test_close(result,'SLE averaged Liouvillian',full(L0),full(H0),...
                  1e-6*norm(H0,1),0,...
                  'the lattice-averaged Liouvillian must be the isotropic lab frame Hamiltonian');

% Check the relaxation superoperator against Redfield theory
result=test_close(result,'SLE relaxation superoperator',full(R_sle),full(R_rf),...
                  1e-4*norm(R_rf,1),0,...
                  'the back-projected superoperator must match Redfield theory to the second order in the anisotropy times correlation time');

end

