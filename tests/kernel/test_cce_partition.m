% Physical grouping, exact group corrections and bounded adaptation tests.
%
%                         result=test_cce_partition()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_partition()

% Build a three-spin conditional problem with an independent exact endpoint
result=new_test_result('kernel/cce_partition','CCE physical partitions',...
                       'Singleton equivalence, exact unions and adaptive held-out gates.');
sys.isotopes={'E','13C','13C','13C'}; sys.magnet=0; sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),diag([0 0 120]),diag([0 0 180]),diag([0 0 270])};
inter.zeeman.units='Hz'; inter.coordinates={[0 0 0],[1 0 0],[2 0 0],[4 0 0]};
inter.coupling.matrix=cell(4);
inter.coupling.matrix{1,2}=[0 0 90;0 0 0;0 0 230];
inter.coupling.matrix{1,3}=[0 0 70;0 0 0;0 0 180];
inter.coupling.matrix{1,4}=[0 0 50;0 0 0;0 0 110];
inter.coupling.matrix{2,3}=30*eye(3); inter.coupling.matrix{3,4}=20*eye(3);
ss=create(sys,inter);
ss.inter.magnet=0; ss.inter.coordinates=inter.coordinates;
ss.inter.zeeman.matrix=cellfun(@(x)2*pi*x,inter.zeeman.matrix,'UniformOutput',false);
ss.inter.coupling.matrix=cell(4);
for j=1:4
    for k=1:4
        if ~isempty(inter.coupling.matrix{j,k}), ss.inter.coupling.matrix{j,k}=2*pi*inter.coupling.matrix{j,k}; end
    end
end
p=parameters(); ordinary=cce(ss,p); full=p; full.order=3; exact=cce(ss,full);

% Singleton partitions have exactly the ordinary CCE subset family
p.expansion='partition'; p.groups={2,3,4}; p.partition=struct('mode','user','max_dimension',8);
singletons=cce(ss,p);
result=test_close(result,'pCCE singleton limit equals ordinary CCE',singletons.value,ordinary.value,0,0,'pCCE singleton limit equals ordinary CCE');
p.groups={[2 3 4]}; p.order=1; whole=cce(ss,p);
result=test_close(result,'One all-bath group gives exact propagation',whole.value,exact.value,2e-12,2e-12,'One all-bath group gives exact propagation');

% Two-group pair corrections retain all cross-group Hamiltonian terms
p.groups={[2 3],4}; p.expansion='disjoint_inter'; paired=cce(ss,p);
result=test_close(result,'Group-order two correction is exact for two groups',paired.value,exact.value,2e-12,2e-12,'Group-order two correction is exact for two groups');
p.expansion='disjoint'; disjoint=cce(ss,p);
result=test_true(result,'Pair correction resolves a nonzero omitted interaction',max(abs(disjoint.value-paired.value))>1e-8,'Pair correction resolves a nonzero omitted interaction');

% Check deterministic coupling and capacity-aware Euclidean grouping
model=cce_prepare(ss,p); graph=cce_graph(model,p); p.partition.mode='coupling';
p.partition.max_dimension=4; a=cce_partition(model,graph,p); b=cce_partition(model,graph,p);
result=test_true(result,'Coupling grouping deterministic with Hilbert ceiling',isequal(a,b)&&all(a.dimensions<=4),'Coupling grouping deterministic with Hilbert ceiling');
p.partition=struct('mode','spatial','max_dimension',4,'group_count',2,'max_rounds',20);
a=cce_partition(model,graph,p);
result=test_true(result,'Spatial grouping preserves a physical partition',a.complete&&isequal(sort([a.groups{:}]),[2 3 4]),'Spatial grouping preserves a physical partition');

% Compare exact pilot patches, freeze groups, then use a distinct held-out state
p=parameters(); p.expansion='disjoint'; model=cce_prepare(ss,p); graph=cce_graph(model,p);
ensemble=cce_states(model,p); held=ensemble.contexts;
held.bath_rho{1}=[1 0;0 0];
p.partition=struct('mode','adaptive','max_dimension',8,'patch_atol',1e-10,...
                   'patch_rtol',1e-10,'max_rounds',5,'max_probes',100,...
                   'patches',{{[2 3 4]}},'pilot_contexts',ensemble.contexts,'heldout_contexts',held);
a=cce_partition(model,graph,p);
result=test_true(result,'Exact adaptive pilot and held-out pass',a.complete&&~isempty(a.diagnostics.heldout_error),'Exact adaptive pilot and held-out pass');
p.partition.max_probes=1; limited=cce_partition(model,graph,p);
result=test_true(result,'Probe ceiling exposes incomplete adaptation',~limited.complete&&~isempty(limited.diagnostics.reason),'Probe ceiling exposes incomplete adaptation');
p.partition.max_probes=100; p.partition.max_dimension=4; limited=cce_partition(model,graph,p);
result=test_true(result,'Exact patch Hilbert ceiling exposes incomplete adaptation',~limited.complete,'Exact patch Hilbert ceiling exposes incomplete adaptation');
fprintf('CCE partition: singleton %.3g, all-group %.3g, pair-group %.3g; adaptive probes %d.\n',...
        max(abs(singletons.value-ordinary.value)),max(abs(whole.value-exact.value)),...
        max(abs(paired.value-exact.value)),a.diagnostics.probes);

end

function p=parameters()

% Declare every mandatory model, state, sequence and numerical choice
p.central_spins=1; p.bath_spins=2:4; p.method='cce'; p.expansion='standard';
p.quantity='coherence'; p.order=2; p.times=linspace(0,0.002,9);
p.central.alpha=[1;0]; p.central.beta=[0;1];
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph.mode='complete'; p.graph.score='coupling_norm'; p.bath_state.mode='mixed_product';
p.meanfield='none'; p.sequence=struct('type','ideal','target',1,'placement','fraction','at',0.5,'axis','x','angle',pi); p.numerics.max_memory_bytes=2^28;
p.numerics.atol=1e-10; p.numerics.rtol=1e-10; p.numerics.reduction='direct'; p.numerics.zero_policy='mask'; p.numerics.propagator='eig';
p.averaging='after_expansion';

end
