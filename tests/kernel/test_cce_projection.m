% Tests central expectations and virtual-level corrections. Syntax:
%
%                        result=test_cce_projection()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_projection()

% State the independent projection targets
result=new_test_result('kernel/cce_projection','CCE central projection',...
                       'Check explicit states, conditional fields and resolved virtual levels.');

% Create a central spin-half with an exactly known longitudinal gap
sys.isotopes={'1H','13C'}; sys.magnet=0; sys.output='hush';
system=create(sys,struct()); gap=1000; coupling=4;
system.inter.zeeman.matrix{1}=gap*eye(3);
system.inter.coupling.matrix{1,2}=diag([coupling 0 0]);
parameters.central_spins=1; parameters.bath_spins=2;
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.method='cce';
parameters.numerics=struct('atol',1e-10,'rtol',1e-12,'max_memory_bytes',1e8);
parameters.central.alpha=[1;0]; parameters.central.beta=[0;1];
model=cce_prepare(system,parameters);
branches=cce_project(model,parameters);
result=test_close(result,'alpha expectation',...
                  branches.alpha.expectations,[0 0 1/2],0,0,'Spin-half up state has positive one-half Sz.');
result=test_close(result,'beta expectation',...
                  branches.beta.expectations,[0 0 -1/2],0,0,'Spin-half down state has negative one-half Sz.');
result=test_true(result,'second order off by default',...
                  ~branches.second_order_enabled&&norm(branches.alpha.susceptibility,'fro')==0,'No implicit perturbative correction.');

% Compare branch Hamiltonians against explicit full-Hilbert projections
[H,meta]=cce_hamiltonian(model,2,struct());
for name={'alpha','beta'}
    branch=branches.(name{1}); embedding=kron(branch.state,eye(2));
    projected=embedding'*H*embedding;
    reference=branch.energy*eye(2)+coupling*branch.expectations(1)*[0 1;1 0]/2;
    result=test_close(result,[name{1} ' conditional projection'],...
                  projected,reference,meta.roundoff_bound,0,'Independent dense product-space projection.');
end

% Enable a central-gap expansion and compare against the exact coupled spectrum
parameters.second_order=struct('enabled',true,'gap_tolerance',1,'overlap_tolerance',1e-12);
branches=cce_project(model,parameters);
result=test_close(result,'virtual coefficient',...
                  branches.alpha.susceptibility(1,1),1/(4*gap),eps,0,'Analytic spin-half transverse susceptibility.');
correction=coupling^2*branches.alpha.susceptibility(1,1)/4;
exact=sort(eig(full(H))); approximate=gap/2+correction;
error_value=abs(approximate-exact(end));
bound=coupling^4/(256*gap^3);
result=test_true(result,'second-order spectrum',...
                  error_value<=bound+eps(gap),'Exact square-root spectrum bounds omitted fourth order.');
fprintf('CCE projection: exact spectral error=%g, fourth-order bound=%g\n',error_value,bound);

% Reject a near-degenerate virtual denominator before any division
parameters.second_order.gap_tolerance=2*gap; rejected=false;
try
    cce_project(model,parameters);
catch
    rejected=true;
end
result=test_true(result,'near-degenerate denominator rejected',...
                  rejected,'An explicit gap guard prevents unresolved perturbation theory.');

% Reject non-eigenstate projectors and gCCE double counting
parameters.second_order.gap_tolerance=1;
parameters.central.alpha=[1;1]/sqrt(2); parameters.central.beta=[1;-1]/sqrt(2);
rejected=false;
try
    cce_project(model,parameters);
catch
    rejected=true;
end
result=test_true(result,'overlap guard',...
                  rejected,'Explicit superposition is not silently reselected.');
parameters.central.alpha=[1;0]; parameters.central.beta=[0;1]; parameters.method='gcce';
rejected=false;
try
    cce_project(model,parameters);
catch
    rejected=true;
end
result=test_true(result,'gCCE correction rejected',...
                  rejected,'Full central dynamics already contains virtual transitions.');

% Verify multi-centre product ordering and a complex orthogonal qubit basis
parameters=rmfield(parameters,'second_order');
parameters.central_spins=[2 1]; parameters.bath_spins=[];
parameters.central.alpha=kron([1;0],[1;1i]/sqrt(2));
parameters.central.beta=kron([0;1],[1;-1i]/sqrt(2));
branches=cce_project(cce_prepare(system,parameters),parameters);
result=test_close(result,'multiple central register',...
                  branches.alpha.expectations,[0 0 1/2;0 1/2 0],8*eps,0,'Complex state expectations preserve declared central order.');

end


