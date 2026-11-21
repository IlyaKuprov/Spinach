% Deterministic serial/parallel scheduling and resource acceptance test. Syntax:
%
%                    result=test_cce_parallel()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
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
sys.isotopes={'13C','13C','13C'};
sys.magnet=0;
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 17]);
inter.coupling.matrix{1,3}=diag([0 0 29]);
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.001 0.004];
p.orientation=[0 0 0];
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.cce=struct('approximation','projected','bath_spins',[2 3],...
             'expansion','standard','order',2,'interaction_cutoff',Inf,...
             'graph',struct('mode','complete','score','coupling_norm'),...
             'meanfield','none','averaging','after_expansion',...
             'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                                 'samples',4,'seed',193),...
             'numerics',struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^26));
[p,~]=cce_options(spin_system,p,@free_evolution,'labframe');

% Measure serial and one-worker execution on the identical sample list
p.cce.parallel=struct('mode','serial','batch_size',3,'max_tasks',16);
started=tic;
[serial,serial_data]=clexpn(spin_system,@free_evolution,p,'labframe');
serial_time=toc(started);
p.cce.parallel.mode='clusters';
started=tic;
[one,one_data]=clexpn(spin_system,@free_evolution,p,'labframe');
one_time=toc(started);
result=test_true(result,'one-worker raw identity',isequaln(serial_data.raw,one_data.raw),...
                 'Indexed raw records are bitwise identical.');
result=test_close(result,'one-worker reduction',one,serial,0,0,...
                  'The final reduction runs only on the client.');

% Close this test's pool before create because create reuses an existing pool
delete(gcp('nocreate'));
sys.parallel={'local',2};
spin_system=create(sys,inter);
pool=gcp('nocreate');
result=test_true(result,'Requested two-worker pool',pool.NumWorkers==2,...
                 'The concurrency and memory checks require two actually running workers.');
started=tic;
[two,two_data]=clexpn(spin_system,@free_evolution,p,'labframe');
two_time=toc(started);
result=test_true(result,'two-worker sample identity',isequal(vertcat(serial_data.ensemble.spin_config.assignment),vertcat(two_data.ensemble.spin_config.assignment)),...
                 'State substreams depend on seed and state ID, never worker rank.');
result=test_true(result,'two-worker raw identity',isequaln(serial_data.raw,two_data.raw),...
                 'Cluster records have a deterministic gather order.');
result=test_close(result,'two-worker reduction',two,serial,0,0,...
                  'The complex result is bitwise identical across worker counts.');
p.cce.parallel.mode='states';
[states,~]=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_close(result,'state scheduling',states,serial,0,0,...
                  'There is one outer parallel layer for either scheduling mode.');

% Compare against independent analytic projected Ising phases
assignments=vertcat(serial_data.ensemble.spin_config.assignment);
magnetic=1.5-assignments;
reference=mean(exp(-1i*(magnetic*[17;29])*2*pi.*p.times),1).';
result=test_close(result,'analytic sampled signal',serial,reference,1e-12,1e-12,...
                  'Projected phases follow the exact diagonal Ising Hamiltonian.');

% Refuse oversized task ledgers and per-worker storage before propagation
p.cce.parallel.max_tasks=1;
[capped,capped_data]=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_true(result,'task cap',isequaln(capped_data.raw,serial_data.raw)&&isequaln(capped,serial),...
                 'The native batch cap limits concurrent work without dropping any physical task.');
p.cce.parallel.max_tasks=16;
p.cce.numerics.max_memory_bytes=2000;
refused=false;
try
    clexpn(spin_system,@free_evolution,p,'labframe');
catch exception
    refused=ismember(exception.identifier,{'cce_clusters:memory',...
                                           'cce_schedule:memory','clexpn:memory'});
end
result=test_true(result,'memory cap',refused,...
                 'Configuration and cluster metadata are bounded before scheduler working matrices.');

% Bound all concurrently active cluster workspaces in the selected batch
p.cce.numerics.max_memory_bytes=2^26;
model=cce_compile(spin_system,p,'labframe');
clusters=serial_data.clusters;
ensemble=serial_data.ensemble;
schedule=struct('pulse_sequence',@free_evolution,'contract',cce_capabilities(@free_evolution));
p.cce.parallel.mode='clusters';
p.cce.parallel.batch_size=2;
count=numel(ensemble.spin_config)*(numel(clusters)+1);
nt=numel(p.times);
dc=prod(model.multiplicities(model.observed_ids));
metadata=whos('model','clusters','ensemble','schedule','p');
resident=3*count*(nt*(25+16*dc^2)+1024)+sum([metadata.bytes]);
dimension=max([clusters.dimension])/dc;
construction=16*(12+3*numel(model.physical_ids))*(dc*dimension)^2;
workspace=construction+16*24*dimension^2;
p.cce.numerics.max_memory_bytes=resident+1.5*workspace;
refused=false;
try
    cce_schedule(model,clusters,ensemble,schedule,p,schedule.contract.evaluator);
catch exception
    if ~strcmp(exception.identifier,'cce_schedule:memory')
        rethrow(exception);
    end
    refused=true;
end
result=test_true(result,'Concurrent worker memory preflight',refused,...
                 'A cap fitting one workspace but not two must reject a two-worker batch before evaluation.');
p.cce.parallel.mode='serial';
accepted=cce_schedule(model,clusters,ensemble,schedule,p,schedule.contract.evaluator);
result=test_true(result,'Single-worker memory preflight',numel(accepted)==numel(ensemble.spin_config),...
                 'The same cap must still admit one serial workspace.');
fprintf('Memory cap %.0f bytes; one workspace %.0f; two-worker refusal %d.\n',...
        p.cce.numerics.max_memory_bytes,workspace,refused);

% Include concurrent Liouville workspaces in the master-equation bound
p.cce.basis.formalism='zeeman-liouv';
p.bath_jumps=struct([]);
p.cce.parallel.mode='clusters';
model=cce_compile(spin_system,p,'labframe');
model.formalism='zeeman-liouv';
metadata=whos('model','clusters','ensemble','schedule','p');
resident=3*count*(nt*(25+16*dc^2)+1024)+sum([metadata.bytes]);
workspace=construction+16*24*dimension^4;
p.cce.numerics.max_memory_bytes=resident+1.5*workspace;
refused=false;
try
    cce_schedule(model,clusters,ensemble,schedule,p,schedule.contract.evaluator);
catch exception
    if ~strcmp(exception.identifier,'cce_schedule:memory')
        rethrow(exception);
    end
    refused=true;
end
result=test_true(result,'Concurrent ME worker memory preflight',refused,...
                 'Each active ME task needs its own squared-Hilbert-space dense workspace.');

% Bound full observed-observed-spin system construction separately from projected propagation
large_sys=struct('isotopes',{repmat({'1H'},1,7)},'magnet',0,'output','hush',...
                 'disable',{{'hygiene'}},'parallel',{{'processes',2}});
large_spin=create(large_sys,struct());
large=p;
large.observed=1:6;
large.cce.bath_spins=7;
large.cce.order=1;
large.times=[0 0.01];
large.cce.basis.formalism='zeeman-hilb';
large=rmfield(large,'bath_jumps');
large.transition=struct('spins',1:6,'alpha',[1;zeros(63,1)],...
                        'beta',[zeros(63,1);1]);
large.cce.bath_state=struct('mode','mixed_product');
large.cce.parallel.batch_size=2;
large.cce.numerics.max_memory_bytes=2^28;
[large,~]=cce_options(large_spin,large,@free_evolution,'labframe');
[~,large_data]=clexpn(large_spin,@free_evolution,large,'labframe');
large_model=cce_compile(large_spin,large,'labframe');
large_clusters=large_data.clusters;
large_ensemble=large_data.ensemble;
for method={'cce','mecce','noise'}

    % Keep an independent high-budget reference for each physical evaluator
    large.cce.approximation='projected';
    large_model.formalism='zeeman-hilb';
    pulse_sequence=@free_evolution;
    if strcmp(method{1},'mecce')
        large.bath_jumps=struct([]);
        large_model.formalism='zeeman-liouv';
    end
    if strcmp(method{1},'noise')
        large=rmfield(large,'bath_jumps');
        large.cce.approximation='generalized';
        pulse_sequence=@noise_correlation;
        large.noise=struct('terms',struct('support',7,'operator',diag([1 -1])/2),...
                           'ordering','ordered','connected',false,'manifold','bath_only');
    end
    large_schedule=struct('pulse_sequence',pulse_sequence,'contract',cce_capabilities(pulse_sequence));
    evaluator=large_schedule.contract.evaluator;
    large_model.parameters=large;
    large.cce.numerics.max_memory_bytes=2^28;
    high=cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large,evaluator);
    metadata=whos('large_model','large_clusters','large_ensemble','large_schedule','large');
    resident=3*2*(2*(25+16*64^2)+1024)+sum([metadata.bytes]);
    if strcmp(method{1},'noise')
        resident=resident+3*2*(32*2+16);
    end
    full_workspace=16*(64*2)^2*(12+3*7);
    if strcmp(method{1},'mecce')

        % Include both retained manifold/spin caches and the active task cache
        empty_cache=16*(3*6*64^2+2);
        full_cache=16*(3*7*128^2+2*2^2);
        resident=resident+3*(empty_cache+full_cache);
        full_workspace=full_workspace+full_cache;

    end
    large.cce.numerics.max_memory_bytes=resident+1.5*full_workspace;
    refused=false;
    try
        cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large,evaluator);
    catch exception
        if ~strcmp(exception.identifier,'cce_schedule:memory')
            rethrow(exception);
        end
        refused=true;
    end
    result=test_true(result,['Full-observed-spin system concurrent workspace ' method{1}],refused,...
                     'Two workers must reserve unprojected observed-observed-spin system storage as well as bath propagation.');

    % The same ceiling must still permit one worker and preserve its numerical result
    large.cce.parallel.mode='serial';
    low=cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large,evaluator);
    result=test_close(result,['Full-observed-spin system serial parity ' method{1}],...
                      low{1}(end).value,high{1}(end).value,0,0,...
                      'The memory guard must reject concurrency, not an individually fitting task.');
    large.cce.parallel.mode='clusters';
    fprintf('CCE full-observed-spin system %s budget %.0f one workspace %.0f refused %d.\n',...
            method{1},large.cce.numerics.max_memory_bytes,full_workspace,refused);

end

% Noise still constructs a full retained Liouville generator before bath restriction
noise_request=large;
noise_request.cce.numerics.max_memory_bytes=2^28;
noise_request.cce.basis.formalism='zeeman-liouv';
noise_request.bath_jumps=struct('support',7,'operator',diag([1 -1])/2,'rate',1);
noise_model=large_model;
noise_model.formalism='zeeman-liouv';
refused=false;
try
    cce_schedule(noise_model,large_clusters,large_ensemble,large_schedule,...
                 noise_request,large_schedule.contract.evaluator);
catch exception
    refused=strcmp(exception.identifier,'cce_schedule:memory');
end
result=test_true(result,'Full Liouville noise memory preflight',refused,...
                 'Bath restriction does not remove the full retained Liouville construction cost.');

% Apply the caller's ceiling to full-observed-spin system projected jump construction too
liouville=large;
liouville.cce.approximation='projected';
liouville=rmfield(liouville,'noise');
liouville.bath_jumps=struct('support',1,'operator',diag([1 -1])/2,'rate',1);
liouville.cce.numerics.max_memory_bytes=2^28;
large_model.formalism='zeeman-liouv';
large_model.parameters=liouville;
[subsystem,~]=cce_subsystem(large_model,7,struct(),[]);
[H,~,~]=cce_generators(subsystem,liouville,[]);
reference=-1i*(kron(eye(2),H.left)-kron(H.right.',eye(2)))+H.dissipator;
liouville.cce.numerics.max_memory_bytes=10000;
refused=false;
try
    cce_generators(subsystem,liouville,[]);
catch exception
    if ~strcmp(exception.identifier,'cce:memory')
        rethrow(exception);
    end
    refused=true;
end
result=test_true(result,'Projected full-observed-spin system Liouvillian guard',refused,...
                 'A projected 4-by-4 generator cannot hide its 128-by-128 physical jump matrices.');
result=test_close(result,'Projected jump independent dephasing',reference,-eye(4)/2,1e-14,1e-14,...
                  'Opposite observed Iz eigenvalues give a bath-independent signal decay rate of one half.');

% Count every observable column in closed and open-system resident records
q=p;
q.cce.bath_state=struct('mode','mixed_product');
q.rho0=eye(2)/2;
q.output='observable';
q.coil=repmat({eye(2)},1,200);
q.times=linspace(0,1,301);
q.cce.parallel.mode='serial';
q.cce.numerics.max_memory_bytes=5000000;
q.bath_jumps=struct([]);
for method={'gcce','megcce'}
    q.cce.approximation='generalized';
    q.cce.basis.formalism='zeeman-hilb';
    if strcmp(method{1},'megcce')
        q.cce.basis.formalism='zeeman-liouv';
    end
    refused=false;
    try
        clexpn(spin_system,@free_evolution,q,'labframe');
    catch exception
        refused=strcmp(exception.identifier,'cce_schedule:memory');
    end
    result=test_true(result,['observable column memory ' method{1}],refused,...
                     'The four time-by-operator records alone exceed five million bytes.');
    fprintf('R7 observable memory %s rejection %d raw minimum %d budget %d\n',...
            method{1},refused,4*25*numel(q.times)*numel(q.coil),q.cce.numerics.max_memory_bytes);
end

% Retained ME manifold, spin, and jump operators accumulate across sampled tasks
for method={'mecce','megcce'}
    retained.observed=1;
    retained.times=[0 1e-3];
    retained.orientation=[0 0 0];
    retained.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
    retained.cce=struct('bath_spins',2,'order',1,'approximation','projected',...
                        'expansion','standard','interaction_cutoff',Inf,...
                        'graph',struct('mode','complete','score','coupling_norm'),...
                        'meanfield','none','averaging','after_expansion',...
                        'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                                            'samples',40,'seed',31),...
                        'numerics',struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',1e8,...
                                          'zero_policy','mask'),...
                        'parallel',struct('mode','serial','max_tasks',80),...
                        'basis',struct('formalism','zeeman-liouv'));
    retained.bath_jumps=repmat(struct('support',2,'operator',[0 1;0 0],'rate',2),1,8);
    if strcmp(method{1},'megcce')
        retained.cce.approximation='generalized';
        retained.output='density';
        retained.rho0=ones(2)/2;
    end
    [retained,~]=cce_options(spin_system,retained,@free_evolution,'labframe');
    [~,retained_data]=clexpn(spin_system,@free_evolution,retained,'labframe');
    retained_model=cce_compile(spin_system,retained,'labframe');
    retained_model.formalism='zeeman-liouv';
    retained_clusters=retained_data.clusters;
    retained_states=retained_data.ensemble;
    retained_schedule=schedule;
    records=retained_data.raw;

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
    components=1;
    if strcmp(method{1},'megcce')
        components=4;
    end
    resident=3*80*(2*(25*components+16*2^2)+1024)+sum([metadata.bytes]);
    dimension=4;
    if strcmp(method{1},'mecce')
        dimension=2;
    end
    workspace=16*4^2*18+16*dimension^4*24;
    retained.cce.numerics.max_memory_bytes=resident+workspace+payload/2;
    refused=false;
    try
        cce_schedule(retained_model,retained_clusters,retained_states,...
                     retained_schedule,retained,retained_schedule.contract.evaluator);
    catch exception
        refused=strcmp(exception.identifier,'cce_schedule:memory');
    end
    result=test_true(result,[method{1} ' retained metadata memory'],refused,...
                     'Every completed task retains its manifold, spin and jump operator arrays.');
    fprintf('%s retained operator payload %.0f bytes; cap %.0f; refused %d.\n',...
            method{1},payload,retained.cce.numerics.max_memory_bytes,refused);
end

% Print the pilot timings without imposing an unjustified speedup target
fprintf('CCE pilot serial %.6g s; one worker %.6g s; two workers %.6g s; discrepancy %.3g.\n',...
        serial_time,one_time,two_time,max(abs(serial-two),[],'all'));
sys.parallel={'local',1};
create(sys,inter);

end


