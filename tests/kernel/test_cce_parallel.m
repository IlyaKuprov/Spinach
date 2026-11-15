% Deterministic serial/parallel scheduling and resource acceptance test.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_parallel()

% Declare the acceptance test and report unavailable parallel execution
result=new_test_result('kernel/cce_parallel','CCE deterministic scheduling',...
                       'Serial, one/two workers, common samples and resource caps.');
if ~license('test','Distrib_Computing_Toolbox')
    error('test_cce_parallel:unavailable',...
          'SKIPPED: Parallel Computing Toolbox unavailable; this harness has no skip status.');
end

% Construct a small analytic Ising bath with deterministic sampled states
sys.isotopes={'13C','13C','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3); inter.coupling.matrix{1,2}=diag([0 0 17]);
inter.coupling.matrix{1,3}=diag([0 0 29]);
spin_system=create(sys,inter);
p=struct('observed_spins',1,'bath_spins',[2 3],'method','cce',...
         'expansion','standard','quantity','signal','times',[0 0.001 0.004],...
         'order',2,'transition',struct('alpha',[1;0],'beta',[0;1]),...
         'assumptions','labframe','orientation',[0 0 0],'interaction_cutoff',Inf,...
         'graph',struct('mode','complete','score','coupling_norm'),...
         'meanfield','none','sequence',[],'averaging','after_expansion',...
         'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                             'samples',4,'seed',193),...
         'numerics',struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^26));

% Measure serial and one-worker execution on the identical sample list
p.parallel=struct('mode','serial','batch_size',3,'max_tasks',16);
started=tic; serial=cce(spin_system,p); serial_time=toc(started);
p.parallel.mode='clusters'; started=tic;
one=cce(spin_system,p); one_time=toc(started);
result=test_true(result,'one-worker raw identity',isequaln(serial.raw,one.raw),...
                 'Indexed raw records are bitwise identical.');
result=test_close(result,'one-worker reduction',one.value,serial.value,0,0,...
                  'The final reduction runs only on the client.');

% Close this test's pool before create because create reuses an existing pool
delete(gcp('nocreate'));
sys.parallel={'local',2}; spin_system=create(sys,inter);
pool=gcp('nocreate');
result=test_true(result,'Requested two-worker pool',pool.NumWorkers==2,...
                 'The concurrency and memory checks require two actually running workers.');
started=tic; two=cce(spin_system,p); two_time=toc(started);
result=test_true(result,'two-worker sample identity',isequal(serial.ensemble.assignments,two.ensemble.assignments),...
                 'State substreams depend on seed and state ID, never worker rank.');
result=test_true(result,'two-worker raw identity',isequaln(serial.raw,two.raw),...
                 'Cluster records have a deterministic gather order.');
result=test_close(result,'two-worker reduction',two.value,serial.value,0,0,...
                  'The complex result is bitwise identical across worker counts.');
p.parallel.mode='states'; states=cce(spin_system,p);
result=test_close(result,'state scheduling',states.value,serial.value,0,0,...
                  'There is one outer parallel layer for either scheduling mode.');

% Compare against independent analytic projected Ising phases
assignments=serial.ensemble.assignments; magnetic=1.5-assignments;
reference=mean(exp(-1i*(magnetic*[17;29])*2*pi.*p.times),1).';
result=test_close(result,'analytic sampled signal',serial.value,reference,1e-12,1e-12,...
                  'Projected phases follow the exact diagonal Ising Hamiltonian.');

% Refuse oversized task ledgers and per-worker storage before propagation
p.parallel.max_tasks=15; refused=false;
try
    cce(spin_system,p); catch exception, refused=strcmp(exception.identifier,'cce_schedule:tasks'); end
result=test_true(result,'task cap',refused,'Sixteen tasks cannot fit a fifteen-task cap.');
p.parallel.max_tasks=16; p.numerics.max_memory_bytes=2000; refused=false;
try
    cce(spin_system,p);
catch exception
    refused=ismember(exception.identifier,{'cce_states:memory','cce_clusters:memory',...
                                          'cce_schedule:memory'});
end
result=test_true(result,'memory cap',refused,...
                 'Configuration and cluster metadata are bounded before scheduler working matrices.');

% Bound all concurrently active cluster workspaces in the selected batch
p.numerics.max_memory_bytes=2^26;
model=cce_prepare(spin_system,p); graph=cce_graph(model,p);
clusters=cce_clusters(model,graph,p); ensemble=cce_states(model,p);
schedule=cce_sequence(model,p); p.parallel.mode='clusters';
p.parallel.batch_size=2; count=numel(ensemble.spin_configs)*(numel(clusters)+1);
nt=numel(schedule.times); dc=prod(model.multiplicities(model.observed_ids));
metadata=whos('model','clusters','ensemble','schedule','p');
resident=3*count*(nt*(25+16*dc^2)+1024)+sum([metadata.bytes]);
dimension=max([clusters.dimension])/dc;
construction=16*(12+3*numel(model.physical_ids)+...
                 2*numel(schedule.events))*(dc*dimension)^2;
workspace=construction+16*(24+numel(schedule.events))*dimension^2;
p.numerics.max_memory_bytes=resident+1.5*workspace; refused=false;
try
    cce_schedule(model,clusters,ensemble,schedule,p,@cce_cluster);
catch exception
    if ~strcmp(exception.identifier,'cce_schedule:memory'), rethrow(exception); end
    refused=true;
end
result=test_true(result,'Concurrent worker memory preflight',refused,...
                 'A cap fitting one workspace but not two must reject a two-worker batch before evaluation.');
p.parallel.mode='serial'; accepted=cce_schedule(model,clusters,ensemble,schedule,p,@cce_cluster);
result=test_true(result,'Single-worker memory preflight',numel(accepted)==numel(ensemble.spin_configs),...
                 'The same cap must still admit one serial workspace.');
fprintf('Memory cap %.0f bytes; one workspace %.0f; two-worker refusal %d.\n',...
        p.numerics.max_memory_bytes,workspace,refused);

% Include concurrent Liouville workspaces in the master-equation bound
p.method='mecce'; p.jumps=struct([]); p.parallel.mode='clusters';
model=cce_prepare(spin_system,p); metadata=whos('model','clusters','ensemble','schedule','p');
resident=3*count*(nt*(25+16*dc^2)+1024)+sum([metadata.bytes]);
workspace=construction+16*(24+numel(schedule.events))*dimension^4;
p.numerics.max_memory_bytes=resident+1.5*workspace; refused=false;
try
    cce_schedule(model,clusters,ensemble,schedule,p,@cce_master);
catch exception
    if ~strcmp(exception.identifier,'cce_schedule:memory'), rethrow(exception); end
    refused=true;
end
result=test_true(result,'Concurrent ME worker memory preflight',refused,...
                 'Each active ME task needs its own squared-Hilbert-space dense workspace.');

% Bound full observed-observed-spin system construction separately from projected propagation
large_sys=struct('isotopes',{repmat({'1H'},1,7)},'magnet',0,'output','hush',...
                 'disable',{{'hygiene'}},'parallel',{{'processes',2}});
large_spin=create(large_sys,struct()); large=p; large.observed_spins=1:6;
large.bath_spins=7; large.order=1; large.times=[0 0.01]; large.method='cce';
large=rmfield(large,'jumps'); large.transition=struct('alpha',[1;zeros(63,1)],...
                                               'beta',[zeros(63,1);1]);
large.bath_state=struct('mode','mixed_product'); large.parallel.batch_size=2;
large.numerics.max_memory_bytes=2^28; large_model=cce_prepare(large_spin,large);
large_graph=cce_graph(large_model,large); large_clusters=cce_clusters(large_model,large_graph,large);
large_ensemble=cce_states(large_model,large); large_schedule=cce_sequence(large_model,large);
for method={'cce','mecce','noise'}

    % Keep an independent high-budget reference for each physical evaluator
    large.method=method{1}; evaluator=@cce_cluster;
    if strcmp(large.method,'mecce'), large.jumps=struct([]); evaluator=@cce_master; end
    if strcmp(large.method,'noise')
        large=rmfield(large,'jumps'); large.quantity='autocorrelation'; evaluator=@cce_noise;
        large.noise=struct('terms',struct('support',7,'operator',diag([1 -1])/2),...
                          'ordering','ordered','connected',false);
    end
    large.numerics.max_memory_bytes=2^28;
    high=cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large,evaluator);
    metadata=whos('large_model','large_clusters','large_ensemble','large_schedule','large');
    resident=3*2*(2*(25+16*64^2)+1024)+sum([metadata.bytes]);
    if strcmp(large.method,'noise'), resident=resident+3*2*(32*2+16); end
    full_workspace=16*(64*2)^2*(12+3*7);
    if strcmp(large.method,'mecce')

        % Include both retained manifold/spin caches and the active task cache
        empty_cache=16*(3*6*64^2+2); full_cache=16*(3*7*128^2+2*2^2);
        resident=resident+3*(empty_cache+full_cache);
        full_workspace=full_workspace+full_cache;
    end
    large.numerics.max_memory_bytes=resident+1.5*full_workspace; refused=false;
    try
        cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large,evaluator);
    catch exception
        if ~strcmp(exception.identifier,'cce_schedule:memory'), rethrow(exception); end
        refused=true;
    end
    result=test_true(result,['Full-observed-spin system concurrent workspace ' large.method],refused,...
                     'Two workers must reserve unprojected observed-observed-spin system storage as well as bath propagation.');

    % The same ceiling must still permit one worker and preserve its numerical result
    large.parallel.mode='serial';
    low=cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large,evaluator);
    result=test_close(result,['Full-observed-spin system serial parity ' large.method],...
                      low{1}(end).value,high{1}(end).value,0,0,...
                      'The memory guard must reject concurrency, not an individually fitting task.');
    large.parallel.mode='clusters';
    fprintf('CCE full-observed-spin system %s budget %.0f one workspace %.0f refused %d.\n',...
            large.method,large.numerics.max_memory_bytes,full_workspace,refused);
end

% Apply the caller's ceiling to full-observed-spin system projected jump construction too
liouville=large; liouville.method='mecce'; liouville.quantity='signal';
liouville=rmfield(liouville,'noise');
liouville.jumps=struct('support',1,'operator',diag([1 -1])/2,'rate',1);
liouville.numerics.max_memory_bytes=2^28;
reference=cce_liouvillian(large_model,7,struct(),liouville);
liouville.numerics.max_memory_bytes=10000; refused=false;
try
    cce_liouvillian(large_model,7,struct(),liouville);
catch exception
    if ~strcmp(exception.identifier,'cce:memory'), rethrow(exception); end
    refused=true;
end
result=test_true(result,'Projected full-observed-spin system Liouvillian guard',refused,...
                 'A projected 4-by-4 generator cannot hide its 128-by-128 physical jump matrices.');
result=test_close(result,'Projected jump independent dephasing',reference,-eye(4)/2,1e-14,1e-14,...
                  'Opposite observed Iz eigenvalues give a bath-independent signal decay rate of one half.');

% Count every observable column in closed and open-system resident records
q=p; q.bath_state=struct('mode','mixed_product'); q.rho0=eye(2)/2;
q.quantity='observable'; q.coil=repmat({eye(2)},1,200);
q.times=linspace(0,1,301); q.parallel.mode='serial';
q.numerics.max_memory_bytes=5000000; q.jumps=struct([]);
for method={'gcce','megcce'}
    q.method=method{1}; refused=false;
    try
        cce(spin_system,q);
    catch exception
        refused=strcmp(exception.identifier,'cce_schedule:memory');
    end
    result=test_true(result,['observable column memory ' method{1}],refused,...
                     'The four time-by-operator records alone exceed five million bytes.');
    fprintf('R7 observable memory %s rejection %d raw minimum %d budget %d\n',...
            method{1},refused,4*25*numel(q.times)*numel(q.coil),q.numerics.max_memory_bytes);
end

% Retained ME manifold, spin and jump operators accumulate across sampled tasks
for method={'mecce','megcce'}
    retained=struct('observed_spins',1,'bath_spins',2,'order',1,...
        'times',[0 1e-3],'method',method{1},'expansion','standard','quantity','signal',...
        'transition',struct('alpha',[1;0],'beta',[0;1]),'rho0',ones(2)/2,...
        'assumptions','labframe','orientation',[0 0 0],'interaction_cutoff',Inf,...
        'graph',struct('mode','complete','score','coupling_norm'),...
        'meanfield','none','sequence',[],'averaging','after_expansion',...
        'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                            'samples',40,'seed',31),...
        'numerics',struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',1e8,...
                          'zero_policy','mask'),...
        'parallel',struct('mode','serial','max_tasks',80));
    retained.jumps=repmat(struct('support',2,'operator',[0 1;0 0],'rate',2),1,8);
    if strcmp(method{1},'megcce'), retained.quantity='density'; end
    retained_model=cce_prepare(spin_system,retained);
    retained_graph=cce_graph(retained_model,retained);
    retained_clusters=cce_clusters(retained_model,retained_graph,retained);
    retained_states=cce_states(retained_model,retained);
    retained_schedule=cce_sequence(retained_model,retained);
    records=cce_schedule(retained_model,retained_clusters,retained_states,...
                         retained_schedule,retained,@cce_master);

    % Bound the actual operator payload without counting MATLAB struct overhead
    payload=0;
    for state=1:numel(records)
        for record=1:numel(records{state})
            meta=records{state}(record).generator;
            elements=numel(meta.h_alpha)+numel(meta.h_beta)+...
                     sum(cellfun(@numel,meta.hamiltonian.operators),'all');
            for jump=meta.jumps
                elements=elements+numel(jump.left)+numel(jump.right);
            end
            payload=payload+16*elements;
        end
    end
    metadata=whos('retained_model','retained_clusters','retained_states',...
                  'retained_schedule','retained');
    components=1; if strcmp(method{1},'megcce'), components=4; end
    resident=3*80*(2*(25*components+16*2^2)+1024)+sum([metadata.bytes]);
    dimension=4; if strcmp(method{1},'mecce'), dimension=2; end
    workspace=16*4^2*18+16*dimension^4*24;
    retained.numerics.max_memory_bytes=resident+workspace+payload/2; refused=false;
    try
        cce_schedule(retained_model,retained_clusters,retained_states,...
                     retained_schedule,retained,@cce_master);
    catch exception
        refused=strcmp(exception.identifier,'cce_schedule:memory');
    end
    result=test_true(result,[method{1} ' retained metadata memory'],refused,...
                     'Every completed task retains its manifold, spin and jump operator arrays.');
    fprintf('%s retained operator payload %.0f bytes; cap %.0f; refused %d.\n',...
            method{1},payload,retained.numerics.max_memory_bytes,refused);
end

% Print the pilot timings without imposing an unjustified speedup target
fprintf('CCE pilot serial %.6g s; one worker %.6g s; two workers %.6g s; discrepancy %.3g.\n',...
        serial_time,one_time,two_time,max(abs(serial.value-two.value),[],'all'));
sys.parallel={'local',1}; create(sys,inter);

end


