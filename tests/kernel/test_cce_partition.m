% Test: kernel/cce_partition | CCE physical partitions
% Physical partitions against an independent dense Kronecker/expm oracle. Syntax:
%
%                    result=test_cce_partition()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_partition()

% Initialise the regression record and a noncommuting three-spin bath
result=new_test_result('kernel/cce_partition','CCE physical partitions',...
                       'Dense exact unions, pair factors and bounded adaptation.');

% Specify noncommuting observed-conditioned tensors in cycles per second
sys.isotopes={'E','13C','13C','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),diag([0 0 120]),diag([0 0 180]),diag([0 0 270])};
inter.zeeman.units='Hz';
inter.coordinates={[0 0 0],[1 0 0],[2 0 0],[4 0 0]};
inter.coupling.matrix=cell(4);
inter.coupling.matrix{1,2}=[0 0 90;0 0 0;0 0 230];
inter.coupling.matrix{1,3}=[0 0 70;0 0 0;0 0 180];
inter.coupling.matrix{1,4}=[0 0 50;0 0 0;0 0 110];
inter.coupling.matrix{2,3}=30*eye(3);
inter.coupling.matrix{3,4}=20*eye(3);
spin_system=create(sys,inter);

% Supply explicit laboratory angular-frequency tensors without field scaling
spin_system.inter.magnet=0;
spin_system.inter.coordinates=inter.coordinates;
spin_system.inter.zeeman.matrix=cellfun(@(x)2*pi*x,inter.zeeman.matrix,'UniformOutput',false);
spin_system.inter.coupling.matrix=cell(4);
for j=1:4
    for k=1:4
        if ~isempty(inter.coupling.matrix{j,k})
            spin_system.inter.coupling.matrix{j,k}=2*pi*inter.coupling.matrix{j,k};
        end
    end
end

parameters=choices();
[parameters,~]=cce_options(spin_system,parameters,@hahn_echo,'labframe');
exact=dense_echo(inter,2:4,parameters.times);
ordinary=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Check the singleton pCCE limit against the identical ordinary CCE family
parameters.cce.expansion='partition';
parameters.cce.groups={2,3,4};
parameters.cce.partition=struct('mode','user','max_dimension',8);
[singletons,singletons_diag]=clexpn(spin_system,@hahn_echo,parameters,'labframe');
result=test_close(result,'Singleton pCCE equals ordinary CCE',...
                  singletons,ordinary,0,0,'Identical subset family.');

% Compare the all-bath group with an independently assembled dense Hamiltonian
parameters.cce.groups={[4 2 3]};
parameters.cce.order=1;
whole=clexpn(spin_system,@hahn_echo,parameters,'labframe');
result=test_close(result,'One group equals dense exact echo',whole,exact,...
                  2e-12,2e-12,'No CCE evaluator is used in the reference.');
result=test_true(result,'Nontrivial dense echo',max(abs(exact-1))>1e-6,...
                 'Reject a commuting or dynamically trivial reference.');

% Check grouped order two and disjoint-inter against exact union propagation
parameters.cce.groups={[2 3],4};
parameters.cce.order=2;
paired=clexpn(spin_system,@hahn_echo,parameters,'labframe');
result=test_close(result,'pCCE group-order two equals dense exact',...
                  paired,exact,2e-12,2e-12,'Both physical groups retained.');
parameters.cce.expansion='disjoint_inter';
parameters.cce.order=1;
inter_result=clexpn(spin_system,@hahn_echo,parameters,'labframe');
result=test_close(result,'Disjoint-inter equals dense exact',...
                  inter_result,exact,2e-12,2e-12,'Pair union keeps cross couplings.');

% Verify disjoint propagation against the product of independent exact groups
parameters.cce.expansion='disjoint';
parameters.cce.order=3;
disjoint=clexpn(spin_system,@hahn_echo,parameters,'labframe');
product=dense_echo(inter,[2 3],parameters.times).*dense_echo(inter,4,parameters.times);
result=test_close(result,'Disjoint equals dense group product',...
                  disjoint,product,2e-12,2e-12,'Disjoint forces group order one.');
result=test_true(result,'Cross-group correction is nonzero',...
                 max(abs(disjoint-exact))>1e-8,'The exact comparison is discriminating.');

% Verify the three-group pair-factor formula independently of the reducer
parameters.cce.groups={2,3,4};
parameters.cce.expansion='disjoint_inter';
inter_result=clexpn(spin_system,@hahn_echo,parameters,'labframe');
pair_product=dense_echo(inter,[2 3],parameters.times).*...
             dense_echo(inter,[2 4],parameters.times).*...
             dense_echo(inter,[3 4],parameters.times)./...
             (dense_echo(inter,2,parameters.times).*...
 dense_echo(inter,3,parameters.times).*...
 dense_echo(inter,4,parameters.times));
result=test_close(result,'Three-group pair-factor identity',...
                  inter_result,pair_product,2e-12,2e-12,...
                  'Each pair contributes F_AB*F0/(F_A*F_B).');

% Check canonical user groups and the singleton default
model=cce_compile(spin_system,parameters,'labframe');
model.pulse_sequence=@hahn_echo;
graph=cce_graph(model,struct('method','cce','graph',parameters.cce.graph));
parameters.cce.groups={[4 2],3};
user=cce_partition(model,graph,parameters);
result=test_true(result,'Canonical physical-ID groups',...
                 isequal(user.groups,{[2 4],3}),'Each user group is sorted, not repartitioned.');
parameters.cce.groups={};
defaults=cce_partition(model,graph,parameters);
result=test_true(result,'Empty groups mean singleton vertices',...
                 isequal(defaults.groups,{2,3,4}),'The documented groups default is honoured.');

% Check strongest coupling merges and deterministic dimension ceilings
parameters.cce.partition.mode='coupling';
parameters.cce.partition.max_dimension=4;
a=cce_partition(model,graph,parameters);
b=cce_partition(model,graph,parameters);
result=test_true(result,'Deterministic strongest-edge grouping',...
                 isequal(a,b)&&isequal(a.groups,{[2 3],4})&&all(a.dimensions<=4),...
                 'The 30 Hz pair precedes the 20 Hz pair.');

% Check capacity-aware spatial assignment and bounded centroid swaps
parameters.cce.partition=struct('mode','spatial','max_dimension',4,...
                                'group_count',2,'max_rounds',10);
a=cce_partition(model,graph,parameters);
b=cce_partition(model,graph,parameters);
result=test_true(result,'Deterministic spatial grouping',...
                 isequal(a,b)&&a.complete&&isequal(a.groups,{[2 3],4}),...
                 'Farthest-first seeds and objective-decreasing swaps.');

% Build distinct fixed pilot and held-out exterior-spin configurations
parameters.cce.expansion='partition';
parameters.cce.order=1;
parameters.cce.groups={2,3,4};
ensemble=singletons_diag.ensemble;
held=ensemble.spin_config;
held.state_id=9001;
held.bath_rho={diag([0 1]),diag([1 0]),diag([0 1])};
parameters.cce.partition=struct('mode','user','adaptive',true,'max_dimension',8,...
                                'patch_atol',1e-10,'patch_rtol',1e-10,...
                                'max_rounds',5,'max_probes',100,'patches',{{2:4}},...
                                'pilot_configs',ensemble.spin_config,'heldout_configs',held);
a=cce_partition(model,graph,parameters);
result=test_true(result,'Adaptive pilot and held-out gates',...
                 a.complete&&isequal(a.groups,{2:4})&&a.diagnostics.rounds>0&&...
                 ~isempty(a.diagnostics.heldout_error),...
                 'Freeze a pilot-selected exact group before held-out testing.');

% Bound diagnostics for both pilot and held-out configuration arrays
for field={'pilot_configs','heldout_configs'}
    q=parameters;
    q.cce.numerics.max_memory_bytes=2^16;
    q.cce.partition.(field{1})=repmat(q.cce.partition.(field{1}),1,4096);
    q.cce.partition.patch_atol=1e6;
    q.cce.partition.max_probes=2;
    refused=false;
    try
        cce_partition(model,graph,q);
    catch exception
        refused=strcmp(exception.identifier,'cce_partition:memory');
    end
    result=test_true(result,['Adaptive ' field{1} ' memory preflight'],refused,...
                     'Pilot and held-out configurations and diagnostics must fit the budget before patch evaluation.');
end

% Restrict every bath-indexed graph matrix to an ordered proper patch
for scoring={'explicit','supplied','resonance'}
    q=parameters;
    q.cce.groups={2:4};
    q.cce.partition.patches={[4 2]};
    q.cce.partition.max_rounds=1;
    q.cce.graph=struct('mode','complete','score','coupling_norm');
    switch scoring{1}
        case 'explicit'
            q.cce.graph.mode='explicit';
            q.cce.graph.adjacency=sparse(~eye(3));
        case 'supplied'
            q.cce.graph.score='supplied';
            q.cce.graph.scores=[0 1 2;1 0 3;2 3 0];
        case 'resonance'
            q.cce.graph.score='resonance';
            q.cce.graph.mismatch=[0 1 2;1 0 3;2 3 0];
            q.cce.graph.mismatch_scale=1;
    end
    passed=false;
    patch_error=Inf;
    try
        probe=cce_partition(model,cce_graph(model,struct('method','cce','graph',q.cce.graph)),q);
        passed=probe.complete;
        patch_error=max(probe.diagnostics.heldout_error);
    catch exception
        fprintf('PATCH_GRAPH %s error=%s\n',scoring{1},exception.identifier);
    end
    result=test_true(result,['Adaptive proper-patch graph ' scoring{1}],passed,...
                     'Full-bath graph matrices follow the declared physical patch ordering.');
    fprintf('PATCH_GRAPH %s passed=%d error=%.3e\n',scoring{1},passed,patch_error);
end

% Expose each exhausted adaptation budget rather than declaring success
parameters.cce.partition.max_probes=1;
limited=cce_partition(model,graph,parameters);
result=test_true(result,'Probe ceiling is incomplete',...
                 ~limited.complete&&limited.diagnostics.probes==1,...
                 'The exhausted probe budget cannot certify the partition.');
parameters.cce.partition.max_probes=100;
parameters.cce.partition.max_dimension=2;
limited=cce_partition(model,graph,parameters);
result=test_true(result,'Exact patch dimension ceiling is incomplete',...
                 ~limited.complete,'An oversized exact patch is not silently dropped.');
parameters.cce.partition.max_dimension=8;
parameters.cce.partition.max_rounds=1;
limited=cce_partition(model,graph,parameters);
result=test_true(result,'Round ceiling is incomplete',...
                 ~limited.complete&&limited.diagnostics.rounds==1,...
                 'One merge cannot make a three-group order-one expansion exact.');

% Require a held-out probe even when the initial pilot already passes
parameters.cce.groups={2:4};
parameters.cce.partition.max_rounds=5;
parameters.cce.partition.max_probes=1;
limited=cce_partition(model,graph,parameters);
result=test_true(result,'Held-out probe cannot be skipped',...
                 ~limited.complete&&contains(limited.diagnostics.reason,'heldout'),...
                 'Pilot success alone is not adaptive completion.');

% Reject an actual held-out discrepancy even when a declared pilot passes
parameters.cce.groups={2,3,4};
parameters.cce.partition.max_probes=100;
parameters.cce.partition.patch_atol=1e-3;
parameters.cce.partition.patch_rtol=0;
held_failure=cce_partition(model,graph,parameters);
result=test_true(result,'Held-out numerical discrepancy is incomplete',...
                 ~held_failure.complete&&contains(held_failure.diagnostics.reason,'heldout'),...
                 'Independent polarized states are not certified by a mixed-state pilot.');
parameters.cce.partition.patch_atol=1e-10;
parameters.cce.partition.patch_rtol=1e-10;
parameters.cce.groups={2:4};

% Reject malformed groups and nonfinite or fractional declared bounds
parameters.cce.partition.adaptive=false;
parameters.cce.partition.max_probes=100;
invalid=parameters;
invalid.cce.groups={[2 3],[3 4]};
refused=false;
try
    cce_partition(model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Overlapping groups rejected',...
                 refused,...
                 'Every physical bath ID occurs exactly once.');
invalid=parameters;
invalid.cce.groups={[2 3]};
refused=false;
try
    cce_partition(model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Missing physical IDs rejected',...
                 refused,...
                 'A partition covers the complete selected bath.');
invalid=parameters;
invalid.cce.partition.max_dimension=NaN;
refused=false;
try
    cce_partition(model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Nonfinite capacity rejected',...
                 refused,...
                 'Declared Hilbert capacities are finite positive integers.');
invalid=parameters;
invalid.cce.partition.adaptive=true;
invalid.cce.partition.max_probes=1.5;
refused=false;
try
    cce_partition(model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Fractional probe ceiling rejected',...
                 refused,...
                 'Probe ceilings count whole exact evaluations.');

% Validate adaptive coordinates, spin_config sets, and tolerances before any probe
invalid=parameters;
invalid.cce.partition.adaptive=true;
invalid.cce.partition.patch_atol=NaN;
refused=false;
try
    cce_partition(model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Nonfinite patch tolerance rejected',...
                 refused,...
                 'NaN tolerances cannot certify a patch.');
invalid=parameters;
invalid.cce.partition.adaptive=true;
invalid.cce.partition.heldout_configs=[];
refused=false;
try
    cce_partition(model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Empty held-out spin_configs rejected',...
                 refused,...
                 'Held-out validation is mandatory for adaptation.');
invalid=parameters;
invalid.cce.groups={};
invalid.cce.partition=struct('mode','spatial','max_dimension',4,...
                             'group_count',2,'max_rounds',1.5);
refused=false;
try
    cce_partition(model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Fractional spatial ceiling rejected',...
                 refused,...
                 'Spatial iteration budgets are whole counts.');
invalid.cce.partition.max_rounds=2;
bad_model=model;
bad_model.coordinates_A(2,1)=NaN;
refused=false;
try
    cce_partition(bad_model,graph,invalid);
catch
    refused=true;
end
result=test_true(result,'Missing spatial coordinate rejected',...
                 refused,...
                 'Spatial grouping requires finite Cartesian coordinates.');

% Keep the empty-bath limit valid for every nonstandard factor family
empty=choices();
empty.cce.bath_spins=[];
empty.cce.order=0;
empty.cce.partition=struct('mode','user','max_dimension',1);
for mode={'partition','disjoint','disjoint_inter'}
    empty.cce.expansion=mode{1};
    empty_result=clexpn(spin_system,@hahn_echo,empty,'labframe');
    result=test_close(result,['Empty bath: ' mode{1}],empty_result,...
                      ones(size(exact)),2e-12,2e-12,'Only the empty-cluster baseline remains.');
end

% Traverse unresolved density reductions until the full three-spin patch is retained
for method={'gcce','megcce'}
    q=choices();
    q.cce.approximation='generalized';
    q.output='density';
    q.rho0=diag([1 0]);
    q.times=[0 1e-4];
    q.cce.expansion='partition';
    q.cce.order=1;
    q.cce.groups={2,3,4};
    if strcmp(method{1},'megcce')
        q.bath_jumps=struct('support',{},'operator',{},'rate',{});
        q.cce.basis.formalism='zeeman-liouv';
        q.cce.numerics=rmfield(q.cce.numerics,'propagator');
    end
    [q,~]=cce_options(spin_system,q,@free_evolution,'labframe');
    local_model=cce_compile(spin_system,q,'labframe');
    local_model.pulse_sequence=@free_evolution;
    if strcmp(method{1},'megcce')
        local_model.formalism='zeeman-liouv';
    end
    states=singletons_diag.ensemble;
    q.cce.partition=struct('mode','user','adaptive',true,'max_dimension',8,...
                           'patch_atol',1e-10,'patch_rtol',1e-10,'max_rounds',2,'max_probes',100,...
                           'patches',{{2:4}},'pilot_configs',states.spin_config,'heldout_configs',states.spin_config);
    local_graph=cce_graph(local_model,struct('method','gcce','graph',q.cce.graph));
    adapted=cce_partition(local_model,local_graph,q);
    result=test_true(result,['Unresolved density intermediate ' method{1}],...
                     adapted.complete&&isequal(adapted.groups,{2:4})&&...
                     adapted.diagnostics.rounds==2&&~isempty(adapted.diagnostics.heldout_error),...
                     'Two bounded merges must reach and validate the exact three-spin density patch.');
    q.cce.partition.max_rounds=1;
    limited=cce_partition(local_model,local_graph,q);
    result=test_true(result,['Unresolved density round ceiling ' method{1}],...
                     ~limited.complete&&numel(limited.groups)==2,...
                     'An unresolved intermediate may advance but cannot be accepted.');
end

% Keep the group partition compatible with generalised observed-density dynamics
general=choices();
general.cce.approximation='generalized';
general.output='density';
general.rho0=ones(2)/2;
general.cce.expansion='partition';
general.cce.order=1;
general.cce.groups={2:4};
general.cce.partition=struct('mode','user','max_dimension',8);
general_result=clexpn(spin_system,@hahn_echo,general,'labframe');
result=test_close(result,'One-group generalized density equals dense echo',...
                  2*general_result(:,3),exact,...
                  2e-12,2e-12,'Generalized CCE retains the complete observed density.');

% Apply exactly the production master-equation evaluator during patch validation
open=choices();
open.cce.basis.formalism='zeeman-liouv';
open.cce.expansion='partition';
open.cce.numerics=rmfield(open.cce.numerics,'propagator');
open.cce.order=1;
open.cce.groups={2:4};
open.cce.partition=parameters.cce.partition;
open.cce.partition.adaptive=true;
open.cce.partition.max_probes=2;
open.bath_jumps=struct('support',2,'operator',diag([1 -1])/2,'rate',100);
[open,~]=cce_options(spin_system,open,@hahn_echo,'labframe');
open_model=model;
open_model.formalism='zeeman-liouv';
open_partition=cce_partition(open_model,graph,open);
result=test_true(result,'Master-equation adaptive evaluator',...
                 open_partition.complete&&open_partition.diagnostics.probes==2,...
                 'Pilot and held-out probes use the same dissipative evaluator.');

% Bypass singular auxiliary factors whenever the complete patch is included
for method={'gcce','megcce'}
    q=choices();
    q.cce.approximation='generalized';
    q.output='density';
    q.rho0=diag([1 0]);
    q.times=[0 1e-4];
    q.cce.expansion='partition';
    q.cce.order=1;
    q.cce.bath_spins=2;
    q.cce.groups={2};
    q.cce.partition=struct('mode','user','max_dimension',2);
    if strcmp(method{1},'megcce')
        q.bath_jumps=struct('support',{},'operator',{},'rate',{});
        q.cce.basis.formalism='zeeman-liouv';
        q.cce.numerics.propagator='stock';
    end
    [reference,reference_diag]=clexpn(spin_system,@free_evolution,q,'labframe');
    [q,~]=cce_options(spin_system,q,@free_evolution,'labframe');
    local_model=cce_compile(spin_system,q,'labframe');
    local_model.pulse_sequence=@free_evolution;
    if strcmp(method{1},'megcce')
        local_model.formalism='zeeman-liouv';
    end
    states=reference_diag.ensemble;
    q.cce.partition.adaptive=true;
    q.cce.partition.patches={2};
    q.cce.partition.patch_atol=1e-10;
    q.cce.partition.patch_rtol=1e-10;
    q.cce.partition.max_rounds=1;
    q.cce.partition.max_probes=10;
    q.cce.partition.pilot_configs=states.spin_config;
    q.cce.partition.heldout_configs=states.spin_config;
    for policy={'mask','error'}
        q.cce.numerics.zero_policy=policy{1};
        accepted=false;
        patch_error=Inf;
        try
            [probe,probe_diag]=clexpn(spin_system,@free_evolution,q,'labframe');
            accepted=all(probe_diag.reduction.valid,'all');
            patch_error=norm(probe(:)-reference(:));
        catch exception
            fprintf('EXACT_PATCH %s %s error=%s\n',method{1},policy{1},exception.identifier);
        end
        result=test_true(result,['Exact adaptive zero baseline ' method{1} ' ' policy{1}],accepted,...
                         'A full-patch density is exact even when auxiliary density entries vanish.');
        result=test_close(result,['Exact adaptive density parity ' method{1} ' ' policy{1}],...
                          patch_error,0,1e-12,1e-12,...
                          'Adaptive validation must agree with the public exact-full-bath bypass.');
        fprintf('EXACT_PATCH %s %s accepted=%d error=%.3e\n',method{1},policy{1},accepted,patch_error);
    end
end

% Reject every incomplete adaptive partition before a one-output production call
for graph_mode={'complete','explicit'}
    for refusal=1:4
        q=choices();
        q.cce.expansion='partition';
        q.cce.order=1;
        q.cce.groups={2,3,4};
        q.cce.graph.mode=graph_mode{1};
        if strcmp(graph_mode{1},'explicit')
            q.cce.graph.adjacency=sparse(~eye(3));
        end
        q.cce.partition=struct('mode','adaptive','max_dimension',8,...
                               'patch_atol',1e-10,'patch_rtol',1e-10,...
                               'max_rounds',5,'max_probes',100,'patches',{{2:4}},...
                               'pilot_configs',ensemble.spin_config,'heldout_configs',held);
        switch refusal
            case 1
                q.cce.partition.max_probes=1;
            case 2
                q.cce.partition.max_dimension=2;
            case 3
                q.cce.partition.max_rounds=1;
            case 4
                q.cce.partition.patch_atol=1e-3;
                q.cce.partition.patch_rtol=0;
        end
        rejected=false;
        try
            clexpn(spin_system,@hahn_echo,q,'labframe');
        catch exception
            rejected=strcmp(exception.identifier,'clexpn:partition');
        end
        result=test_true(result,['Incomplete production partition ' graph_mode{1} ' ' num2str(refusal)],...
                         rejected,'Pilot, held-out, capacity, and iteration failures cannot become ordinary answers.');
    end
end

% Print quantitative evidence without adding nonstandard result fields
fprintf('PARTITION_METRICS dense_error=%.3e pair_error=%.3e disjoint_error=%.3e cross_effect=%.3e adaptive_rounds=%d probes=%d heldout=%.3e\n',...
        max(abs(whole-exact)),max(abs(paired-exact)),...
        max(abs(disjoint-product)),max(abs(disjoint-exact)),...
        a.diagnostics.rounds,a.diagnostics.probes,max(a.diagnostics.heldout_error));
fprintf('PARTITION_HELDOUT pilot_error=%.3e heldout_error=%.3e tolerance=%.3e\n',...
        max(held_failure.diagnostics.pilot_error),...
        max(held_failure.diagnostics.heldout_error),1e-3);

end

% Specify the physical and numerical test parameters
function parameters=choices()

% Declare physical, state, sequence, and numerical choices explicitly
parameters.observed=1;
parameters.cce.bath_spins=2:4;
parameters.cce.approximation='projected';
parameters.cce.expansion='standard';
parameters.output='signal';
parameters.cce.order=2;
parameters.times=linspace(0,0.002,9);
parameters.transition.spins=1;
parameters.transition.alpha=[1;0];
parameters.transition.beta=[0;1];
parameters.orientation=[0 0 0];
parameters.cce.interaction_cutoff=Inf;
parameters.cce.graph.mode='complete';
parameters.cce.graph.score='coupling_norm';
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.meanfield='none';
parameters.cce.numerics=struct('max_memory_bytes',2^28,'atol',1e-10,'rtol',1e-10,...
                               'zero_policy','mask','propagator','eig');
parameters.cce.averaging='after_expansion';

end

% Evaluate an independent dense Hahn-echo reference
function value=dense_echo(inter,sites,times)

% Assemble independent Cartesian spin-half operators in full cluster Hilbert space
ids=[1 sites];
count=numel(ids);
dim=2^count;
spin_ops={[0 1;1 0]/2,[0 -1i;1i 0]/2,diag([1 -1])/2};
ops=cell(count,3);
for n=1:count
    for axis=1:3
        ops{n,axis}=kron(kron(eye(2^(n-1)),spin_ops{axis}),eye(2^(count-n)));
    end
end

% Add the specified local fields and every physical pair exactly once
hamiltonian=zeros(dim);
for n=1:count
    hamiltonian=hamiltonian+2*pi*inter.zeeman.matrix{ids(n)}(3,3)*ops{n,3};
    for k=n+1:count
        tensor=inter.coupling.matrix{ids(n),ids(k)}.';
        if isempty(tensor)
            continue;
        end
        for axis=1:3
            for other=1:3
                hamiltonian=hamiltonian+2*pi*tensor(axis,other)*ops{n,axis}*ops{k,other};
            end
        end
    end
end

% Propagate an observed superposition and maximally mixed bath through a Hahn echo
db=dim/2;
density=kron(ones(2)/2,eye(db)/db);
pulse=expm(-1i*pi*ops{1,1});
value=zeros(numel(times),1);
for k=1:numel(times)
    free=expm(-1i*hamiltonian*times(k)/2);
    propagator=free*pulse*free;
    final=propagator*density*propagator';
    value(k)=2*trace(final(1:db,db+1:end));
end

end


