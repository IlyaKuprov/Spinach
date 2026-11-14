% Physical partitions against an independent dense Kronecker/expm oracle.
% Syntax:
%
%                         result=test_cce_partition()
%
% Outputs:
%     result - standard Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_partition()

% Initialise the regression record and a noncommuting three-spin bath
result=new_test_result('kernel/cce_partition','CCE physical partitions',...
                       'Dense exact unions, pair factors and bounded adaptation.');
[spin_system,inter]=fixture(); parameters=choices();
exact=dense_echo(inter,2:4,parameters.times);
ordinary=cce(spin_system,parameters);

% Check the singleton pCCE limit against the identical ordinary CCE family
parameters.expansion='partition'; parameters.groups={2,3,4};
parameters.partition=struct('mode','user','max_dimension',8);
singletons=cce(spin_system,parameters);
result=test_close(result,'Singleton pCCE equals ordinary CCE',...
                  singletons.value,ordinary.value,0,0,'Identical subset family.');

% Compare the all-bath group with an independently assembled dense Hamiltonian
parameters.groups={[4 2 3]}; parameters.order=1;
whole=cce(spin_system,parameters);
result=test_close(result,'One group equals dense exact echo',whole.value,exact,...
                  2e-12,2e-12,'No CCE evaluator is used in the reference.');
result=test_true(result,'Nontrivial dense echo',max(abs(exact-1))>1e-6,...
                 'Reject a commuting or dynamically trivial reference.');

% Check grouped order two and disjoint-inter against exact union propagation
parameters.groups={[2 3],4}; parameters.order=2;
paired=cce(spin_system,parameters);
result=test_close(result,'pCCE group-order two equals dense exact',...
                  paired.value,exact,2e-12,2e-12,'Both physical groups retained.');
parameters.expansion='disjoint_inter'; parameters.order=1;
inter_result=cce(spin_system,parameters);
result=test_close(result,'Disjoint-inter equals dense exact',...
                  inter_result.value,exact,2e-12,2e-12,'Pair union keeps cross couplings.');

% Verify disjoint propagation against the product of independent exact groups
parameters.expansion='disjoint'; parameters.order=3;
disjoint=cce(spin_system,parameters);
product=dense_echo(inter,[2 3],parameters.times).*dense_echo(inter,4,parameters.times);
result=test_close(result,'Disjoint equals dense group product',...
                  disjoint.value,product,2e-12,2e-12,'Disjoint forces group order one.');
result=test_true(result,'Cross-group correction is nonzero',...
                 max(abs(disjoint.value-exact))>1e-8,'The exact comparison is discriminating.');

% Verify the three-group pair-factor formula independently of the reducer
parameters.groups={2,3,4}; parameters.expansion='disjoint_inter';
inter_result=cce(spin_system,parameters);
pair_product=dense_echo(inter,[2 3],parameters.times).*...
             dense_echo(inter,[2 4],parameters.times).*...
             dense_echo(inter,[3 4],parameters.times)./...
             (dense_echo(inter,2,parameters.times).*...
              dense_echo(inter,3,parameters.times).*...
              dense_echo(inter,4,parameters.times));
result=test_close(result,'Three-group pair-factor identity',...
                  inter_result.value,pair_product,2e-12,2e-12,...
                  'Each pair contributes F_AB*F0/(F_A*F_B).');

% Check canonical user groups and the singleton default
model=cce_prepare(spin_system,parameters); graph=cce_graph(model,parameters);
parameters.groups={[4 2],3}; user=cce_partition(model,graph,parameters);
result=test_true(result,'Canonical physical-ID groups',...
                 isequal(user.groups,{[2 4],3}),'Each user group is sorted, not repartitioned.');
parameters.groups={}; defaults=cce_partition(model,graph,parameters);
result=test_true(result,'Empty groups mean singleton vertices',...
                 isequal(defaults.groups,{2,3,4}),'The documented groups default is honoured.');

% Check strongest coupling merges and deterministic dimension ceilings
parameters.partition.mode='coupling'; parameters.partition.max_dimension=4;
a=cce_partition(model,graph,parameters); b=cce_partition(model,graph,parameters);
result=test_true(result,'Deterministic strongest-edge grouping',...
                 isequal(a,b)&&isequal(a.groups,{[2 3],4})&&all(a.dimensions<=4),...
                 'The 30 Hz pair precedes the 20 Hz pair.');

% Check capacity-aware spatial assignment and bounded centroid swaps
parameters.partition=struct('mode','spatial','max_dimension',4,...
                            'group_count',2,'max_rounds',10);
a=cce_partition(model,graph,parameters); b=cce_partition(model,graph,parameters);
result=test_true(result,'Deterministic spatial grouping',...
                 isequal(a,b)&&a.complete&&isequal(a.groups,{[2 3],4}),...
                 'Farthest-first seeds and objective-decreasing swaps.');

% Build distinct fixed pilot and held-out exterior-state contexts
parameters.expansion='partition'; parameters.order=1; parameters.groups={2,3,4};
ensemble=cce_states(model,parameters); held=ensemble.contexts;
held.id=uint64(9001); held.external_m=[-0.5 0.5 -0.5];
parameters.partition=struct('mode','user','adaptive',true,'max_dimension',8,...
                            'patch_atol',1e-10,'patch_rtol',1e-10,...
                            'max_rounds',5,'max_probes',100,'patches',{{2:4}},...
                            'pilot_contexts',ensemble.contexts,'heldout_contexts',held);
a=cce_partition(model,graph,parameters);
result=test_true(result,'Adaptive pilot and held-out gates',...
                 a.complete&&isequal(a.groups,{2:4})&&a.diagnostics.rounds>0&&...
                 ~isempty(a.diagnostics.heldout_error),...
                 'Freeze a pilot-selected exact group before held-out testing.');

% Expose each exhausted adaptation budget rather than declaring success
parameters.partition.max_probes=1; limited=cce_partition(model,graph,parameters);
result=test_true(result,'Probe ceiling is incomplete',...
                 ~limited.complete&&limited.diagnostics.probes==1,...
                 'The exhausted probe budget cannot certify the partition.');
parameters.partition.max_probes=100; parameters.partition.max_dimension=2;
limited=cce_partition(model,graph,parameters);
result=test_true(result,'Exact patch dimension ceiling is incomplete',...
                 ~limited.complete,'An oversized exact patch is not silently dropped.');
parameters.partition.max_dimension=8; parameters.partition.max_rounds=1;
limited=cce_partition(model,graph,parameters);
result=test_true(result,'Round ceiling is incomplete',...
                 ~limited.complete&&limited.diagnostics.rounds==1,...
                 'One merge cannot make a three-group order-one expansion exact.');

% Require a held-out probe even when the initial pilot already passes
parameters.groups={2:4}; parameters.partition.max_rounds=5;
parameters.partition.max_probes=1; limited=cce_partition(model,graph,parameters);
result=test_true(result,'Held-out probe cannot be skipped',...
                 ~limited.complete&&contains(limited.diagnostics.reason,'heldout'),...
                 'Pilot success alone is not adaptive completion.');

% Reject malformed groups and nonfinite or fractional declared bounds
parameters.partition.adaptive=false; parameters.partition.max_probes=100;
invalid=parameters; invalid.groups={[2 3],[3 4]};
result=test_true(result,'Overlapping groups rejected',...
                 rejected(@()cce_partition(model,graph,invalid)),...
                 'Every physical bath ID occurs exactly once.');
invalid=parameters; invalid.groups={[2 3]};
result=test_true(result,'Missing physical IDs rejected',...
                 rejected(@()cce_partition(model,graph,invalid)),...
                 'A partition covers the complete selected bath.');
invalid=parameters; invalid.partition.max_dimension=NaN;
result=test_true(result,'Nonfinite capacity rejected',...
                 rejected(@()cce_partition(model,graph,invalid)),...
                 'Declared Hilbert capacities are finite positive integers.');
invalid=parameters; invalid.partition.adaptive=true;
invalid.partition.max_probes=1.5;
result=test_true(result,'Fractional probe ceiling rejected',...
                 rejected(@()cce_partition(model,graph,invalid)),...
                 'Probe ceilings count whole exact evaluations.');

% Print quantitative evidence without adding nonstandard result fields
fprintf('PARTITION_METRICS dense_error=%.3e pair_error=%.3e disjoint_error=%.3e cross_effect=%.3e adaptive_rounds=%d probes=%d heldout=%.3e\n',...
        max(abs(whole.value-exact)),max(abs(paired.value-exact)),...
        max(abs(disjoint.value-product)),max(abs(disjoint.value-exact)),...
        a.diagnostics.rounds,a.diagnostics.probes,max(a.diagnostics.heldout_error));

end

function [spin_system,inter]=fixture()

% Specify noncommuting central-conditioned tensors in cycles per second
sys.isotopes={'E','13C','13C','13C'}; sys.magnet=0;
sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),diag([0 0 120]),diag([0 0 180]),diag([0 0 270])};
inter.zeeman.units='Hz'; inter.coordinates={[0 0 0],[1 0 0],[2 0 0],[4 0 0]};
inter.coupling.matrix=cell(4);
inter.coupling.matrix{1,2}=[0 0 90;0 0 0;0 0 230];
inter.coupling.matrix{1,3}=[0 0 70;0 0 0;0 0 180];
inter.coupling.matrix{1,4}=[0 0 50;0 0 0;0 0 110];
inter.coupling.matrix{2,3}=30*eye(3); inter.coupling.matrix{3,4}=20*eye(3);
spin_system=create(sys,inter);

% Supply explicit laboratory angular-frequency tensors without field scaling
spin_system.inter.magnet=0; spin_system.inter.coordinates=inter.coordinates;
spin_system.inter.zeeman.matrix=cellfun(@(x)2*pi*x,inter.zeeman.matrix,'UniformOutput',false);
spin_system.inter.coupling.matrix=cell(4);
for j=1:4
    for k=1:4
        if ~isempty(inter.coupling.matrix{j,k})
            spin_system.inter.coupling.matrix{j,k}=2*pi*inter.coupling.matrix{j,k};
        end
    end
end

end

function parameters=choices()

% Declare physical, state, sequence, and numerical choices explicitly
parameters.central_spins=1; parameters.bath_spins=2:4;
parameters.method='cce'; parameters.expansion='standard';
parameters.quantity='coherence'; parameters.order=2;
parameters.times=linspace(0,0.002,9);
parameters.central.alpha=[1;0]; parameters.central.beta=[0;1];
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf; parameters.graph.mode='complete';
parameters.graph.score='coupling_norm'; parameters.bath_state.mode='mixed_product';
parameters.meanfield='none';
parameters.sequence=struct('type','ideal','target',1,'placement','fraction',...
                           'at',0.5,'axis','x','angle',pi);
parameters.numerics=struct('max_memory_bytes',2^28,'atol',1e-10,'rtol',1e-10,...
                           'reduction','direct','zero_policy','mask','propagator','eig');
parameters.averaging='after_expansion';

end

function value=dense_echo(inter,sites,times)

% Assemble independent Cartesian spin-half operators in full cluster Hilbert space
ids=[1 sites]; count=numel(ids); dim=2^count;
pauli={[0 1;1 0]/2,[0 -1i;1i 0]/2,diag([1 -1])/2};
ops=cell(count,3);
for n=1:count
    for axis=1:3
        ops{n,axis}=kron(kron(eye(2^(n-1)),pauli{axis}),eye(2^(count-n)));
    end
end

% Add the specified local fields and every physical pair exactly once
hamiltonian=zeros(dim);
for n=1:count
    hamiltonian=hamiltonian+2*pi*inter.zeeman.matrix{ids(n)}(3,3)*ops{n,3};
    for k=n+1:count
        tensor=inter.coupling.matrix{ids(n),ids(k)}.';
        if isempty(tensor), continue; end
        for axis=1:3
            for other=1:3
                hamiltonian=hamiltonian+2*pi*tensor(axis,other)*ops{n,axis}*ops{k,other};
            end
        end
    end
end

% Propagate a central superposition and maximally mixed bath through a Hahn echo
db=dim/2; density=kron(ones(2)/2,eye(db)/db);
pulse=expm(-1i*pi*ops{1,1}); value=zeros(numel(times),1);
for k=1:numel(times)
    free=expm(-1i*hamiltonian*times(k)/2); propagator=free*pulse*free;
    final=propagator*density*propagator'; value(k)=2*trace(final(1:db,db+1:end));
end

end

function outcome=rejected(operation)

% Accept only a genuinely thrown validation error
outcome=false;
try
    operation();
catch
    outcome=true;
end

end

