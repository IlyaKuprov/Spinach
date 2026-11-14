% Density physicality and singular-baseline exact bypass tests. Syntax:
%
%                         result=test_gcce_density()
%
% ilya.kuprov@weizmann.ac.il

function result=test_gcce_density()

% State the no-repair requirement for elementwise density expansions
result=new_test_result('kernel/gcce_density','gCCE density diagnostics',...
                       'Baseline zeros mask truncated ratios and full order bypasses them.');
sys.isotopes={'1H','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3)}; inter.coupling.matrix=cell(2);
spin_system=create(sys,inter); spin_system.inter.zeeman.matrix{1}(3,1)=100;
spin_system.inter.coupling.matrix{1,2}=diag([31 43 57]);
p.central_spins=1; p.bath_spins=2; p.method='gcce'; p.expansion='standard';
p.quantity='density'; p.central.rho=[1 0;0 0]; p.assumptions='labframe';
p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm');
p.meanfield='none'; p.times=[0 pi/200 pi/100]; p.sequence=[];
p.bath_state=struct('mode','mixed_product');
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

% Exact full-bath density must stay finite at zero central-only elements
p.order=1; actual=cce(spin_system,p);
result=test_true(result,'full-order bypass finite',all(actual.valid,'all'),...
                 'Exact reduced density does not require a baseline division.');
result=test_close(result,'initial density',actual.density(:,:,1),p.central.rho,1e-12,0,...
                  'Zero duration returns the supplied central density.');
result=test_true(result,'unit trace',max(actual.diagnostics.trace_error)<1e-12,...
                 'Unitary full-system propagation preserves partial-trace normalization.');
result=test_true(result,'Hermitian',max(actual.diagnostics.hermiticity_error)<1e-12,...
                 'No Hermitian repair is applied.');
result=test_true(result,'positive',min(actual.diagnostics.minimum_eigenvalue)>-1e-12,...
                 'Exact partial trace is positive without eigenvalue clipping.');

% Force a truncated family with a zero baseline and retain the raw record
p.bath_spins=2; p.order=0; truncated=cce(spin_system,p);
result=test_true(result,'central baseline preserved',all(isfinite(truncated.value),'all'),...
                 'Order zero returns the direct central-only density.');
p.order=1; p.bath_spins=2; p.central_spins=1;
model=cce_prepare(spin_system,p); graph=cce_graph(model,p); clusters=cce_clusters(model,graph,p);
records=actual.raw{1}; reduced=cce_reduce(records,clusters,p);
result=test_true(result,'elementwise baseline zeros masked',any(~reduced.valid,'all'),...
                 'The reducer alone must not replace undefined ratios by one.');
result=test_true(result,'raw density retained',all(isfinite(records(2).value),'all'),...
                 'Masking factors does not destroy direct raw density information.');

% Reject nonpositive densities and unresolved initial scalar normalization
bad=p; bad.central.rho=diag([1.1 -0.1]); caught=false;
try
    cce_states(model,bad);
catch
    caught=true;
end
result=test_true(result,'negative input rejected',caught,'State positivity is checked.');
bad=p; bad.quantity='coherence'; bad.normalized=true;
bad.central.alpha=[1;0]; bad.central.beta=[0;1]; caught=false;
try
    cce_observe(bad.central.rho,bad);
catch
    caught=true;
end
result=test_true(result,'zero normalization rejected',caught,...
                 'Initial zero coherence has no normalized scalar interpretation.');

% Reproduce fixed and sampled product states without changing the global stream
sample=p; sample.bath_state=struct('mode','sampled','samples',5,'seed',17,...
                                  'basis','uniform','probabilities',{{[1;0]}});
sample.averaging='after_expansion'; before=rng; ensemble=cce_states(model,sample);
after=rng; again=cce_states(model,sample);
result=test_true(result,'global RNG unchanged',isequal(before,after),...
                 'Sampling owns counter-based Threefry substreams.');
result=test_true(result,'fixed probability state',all(ensemble.assignments==1,'all'),...
                 'Probability-one states are not resampled into forbidden levels.');
result=test_true(result,'reproducible assignments',isequal(ensemble.assignments,again.assignments),...
                 'Explicit state IDs determine the stream independently of workers.');
sample.bath_state.probabilities='uniform'; sample.meanfield='static';
ensemble=cce_states(model,sample); schedule=cce_sequence(model,sample);
small=cce_cluster(model,struct('id',0,'sites',[]),ensemble.contexts(1),schedule,sample);
context=ensemble.contexts(1); field=model.retained{1,2}*context.expectations.';
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2;
h=100*x+field(1)*x+field(2)*y+field(3)*z;
reference=zeros(numel(p.times),4);
for k=1:numel(p.times)
    u=expm(-1i*h*p.times(k)); rho=u*p.central.rho*u'; reference(k,:)=rho(:).';
end
result=test_close(result,'static exterior exact',small.value,reference,1e-12,1e-12,...
                  'Independent central-only Hamiltonian with the same exterior assignment.');

% Exercise both explicitly named sampling estimators at exact full order
sample.quantity='density'; sample.order=1;
sample.averaging='after_expansion'; after=cce(spin_system,sample);
sample.averaging='internal_then_expansion'; internal=cce(spin_system,sample);
result=test_close(result,'full-order sampling estimators',after.value,internal.value,...
                  1e-12,1e-12,'Exact raw density is linear in the initial state.');
result=test_true(result,'sample errors finite',all(isfinite(after.abs_error_bound),'all'),...
                 'Weighted sums retain explicit forward bounds.');

% Distinguish selective qubit pulses from full spin-one rotations
toy=model; toy.physical_ids=1; toy.central_ids=1; toy.multiplicities=3;
sequence=struct('type','ideal','target',1,'placement','fraction','at',0.5,...
                'axis','x','angle',pi); sample.sequence=sequence;
compiled=cce_sequence(toy,sample);
result=test_close(result,'explicit pulse timestamps',compiled.positions,p.times/2,0,0,...
                  'Public time is total duration rather than half echo time.');
a=[1;0;0]; b=[0;1;0]; sx=[0 1 0;1 0 1;0 1 0]/sqrt(2);
result=test_true(result,'spin-one selective distinction',...
    norm(expm(-1i*pi*sx)-expm(-1i*pi*(a*b'+b*a')/2),'fro')>1,...
    'A full physical spin-one pi rotation is not a selective two-state gate.');

result.metrics.trace_error=max(actual.diagnostics.trace_error);
result.metrics.hermiticity_error=max(actual.diagnostics.hermiticity_error);
result.metrics.minimum_eigenvalue=min(actual.diagnostics.minimum_eigenvalue);
result.metrics.masked_elements=nnz(~reduced.valid);

end

