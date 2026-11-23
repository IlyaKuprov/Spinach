% Test: kernel/cce_parallel | CCE deterministic scheduling
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
sys.parallel={'processes',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3);
inter.coupling.matrix{1,2}=diag([0 0 17]);
inter.coupling.matrix{1,3}=diag([0 0 29]);
delete(gcp('nocreate'));
sys.tols.cce_atol=1e-8;
sys.tols.cce_rtol=1e-8;
sys.tols.cce_mem=2^26;
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.001 0.004];
p.orientation=[0 0 0];
p.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
p.cce=struct('approximation','projected','bath_spins',[2 3],...
             'expansion','standard','order',2,'interaction_cutoff',Inf,...
             'meanfield','none','averaging','after_expansion');
p.bath_samples=4;
p.bath_seed=193;
p=cce_options(spin_system,p,@free_evolution,'labframe');

% Measure serial and one-worker execution on the identical sample list
delete(gcp('nocreate'));
started=tic;
[serial,serial_data]=clexpn(spin_system,@free_evolution,p,'labframe');
serial_time=toc(started);
create(sys,inter);
started=tic;
[one,one_data]=clexpn(spin_system,@free_evolution,p,'labframe');
one_time=toc(started);
result=test_true(result,'one-worker raw identity',isequaln(serial_data.raw,one_data.raw),...
                 'Indexed raw records are bitwise identical.');
result=test_close(result,'one-worker reduction',one,serial,0,0,...
                  'The final reduction runs only on the client.');

% Close this test's pool before create because create reuses an existing pool
delete(gcp('nocreate'));
sys.parallel={'processes',2};
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

% Compare against independent analytic projected Ising phases
assignments=vertcat(serial_data.ensemble.spin_config.assignment);
magnetic=1.5-assignments;
reference=mean(exp(-1i*(magnetic*[17;29])*2*pi.*p.times),1).';
result=test_close(result,'analytic sampled signal',serial,reference,1e-12,1e-12,...
                  'Projected phases follow the exact diagonal Ising Hamiltonian.');

% Refuse a cap smaller than the ledger and a single task workspace
tight=spin_system;
tight.tols.cce_mem=2000;
refused=false;
try
    clexpn(tight,@free_evolution,p,'labframe');
catch exception
    refused=ismember(exception.identifier,{'cce_clusters:memory',...
                                           'cce_schedule:memory','clexpn:memory'});
end
result=test_true(result,'memory cap',refused,...
                 'Configuration and cluster metadata are bounded before scheduler working matrices.');

% Fit only one workspace on a two-worker pool, in each formalism
clusters=serial_data.clusters;
ensemble=serial_data.ensemble;
schedule=@free_evolution;
count=numel(ensemble.spin_config)*(numel(clusters)+1);
nt=numel(p.times);
dc=prod(spin_system.comp.mults(p.observed));
for formalism={'zeeman-hilb','zeeman-liouv'}
    source=spin_system;
    source.bas.formalism=formalism{1};
    model_p=cce_options(source,p,@cce_compile,'labframe');
    model=cce_compile(source,model_p,'labframe');
    high=cce_schedule(model,clusters,ensemble,schedule,p);
    metadata=whos('model','clusters','ensemble','p');
    resident=3*count*(nt*(25+16*dc^2)+1024)+sum([metadata.bytes]);
    bath_dim=max([clusters.dimension])/dc;
    nspins=numel(p.observed)+max(arrayfun(@(c)numel(c.sites),clusters));
    workspace=16*(12+3*nspins)*(dc*bath_dim)^2+16*24*bath_dim^2;
    if strcmp(formalism{1},'zeeman-liouv')
        empty_cache=16*(3*dc^2+2)+16;
        full_cache=16*(3*nspins*(dc*bath_dim)^2+2*bath_dim^2)+16*nspins;
        resident=resident+3*numel(ensemble.spin_config)*(empty_cache+numel(clusters)*full_cache);
        workspace=workspace+full_cache+16*24*(bath_dim^4-bath_dim^2);
    end
    model.tols.cce_mem=resident+1.5*workspace;
    accepted=cce_schedule(model,clusters,ensemble,schedule,p);
    result=test_true(result,['Single-workspace parity ' formalism{1}],...
                     isequaln(accepted,high),...
                     'The cap admits one workspace without changing the raw records.');
end

% Bound full observed-observed-spin system construction separately from projected propagation
large_sys=struct('isotopes',{repmat({'1H'},1,7)},'magnet',0,'output','hush',...
                 'disable',{{'hygiene'}},'parallel',{{'processes',2}});
large_sys.tols=sys.tols;
large_sys.tols.cce_mem=2^28;
large_spin=create(large_sys,struct());
large=rmfield(p,{'bath_samples','bath_seed'});
large.observed=1:6;
large.cce.bath_spins=7;
large.cce.order=1;
large.times=[0 0.01];
large.transition=struct('spins',1:6,'alpha',[1;zeros(63,1)],...
                        'beta',[zeros(63,1);1]);
large=cce_options(large_spin,large,@free_evolution,'labframe');
[~,large_data]=clexpn(large_spin,@free_evolution,large,'labframe');
large_model_p=cce_options(large_spin,large,@cce_compile,'labframe');
large_model=cce_compile(large_spin,large_model_p,'labframe');
large_clusters=large_data.clusters;
large_ensemble=large_data.ensemble;
for method={'cce','mecce','noise'}

    % Keep an independent high-budget reference for each physical evaluator
    large.cce.approximation='projected';
    large_model.bas.formalism='zeeman-hilb';
    pulse_sequence=@free_evolution;
    if strcmp(method{1},'mecce')
        large.bath_jumps=struct([]);
        large_model.bas.formalism='zeeman-liouv';
    end
    if strcmp(method{1},'noise')
        large=rmfield(large,'bath_jumps');
        large.cce.approximation='generalized';
        pulse_sequence=@noise_correlation;
        large.noise=struct('terms',struct('support',7,'operator',diag([1 -1])/2),...
                           'ordering','ordered','connected',false,'manifold','bath_only');
    end
    large_schedule=pulse_sequence;
    large_model.tols.cce_mem=2^28;
    high=cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large);
    metadata=whos('large_model','large_clusters','large_ensemble','large');
    resident=3*2*(2*(25+16*64^2)+1024)+sum([metadata.bytes]);
    if strcmp(method{1},'noise')
        resident=resident+3*2*(32*2+16);
    end
    full_workspace=16*(64*2)^2*(12+3*7)+16*24*2^2;
    if strcmp(method{1},'mecce')

        % Include both retained manifold/spin caches and the active task cache
        empty_cache=16*(3*6*64^2+2)+16*6;
        full_cache=16*(3*7*128^2+2*2^2)+16*7;
        resident=resident+3*(empty_cache+full_cache);
        full_workspace=full_workspace+full_cache+16*24*(2^4-2^2);

    end
    large_model.tols.cce_mem=resident+0.5*full_workspace;
    refused=false;
    try
        cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large);
    catch exception
        if ~strcmp(exception.identifier,'cce_schedule:memory')
            rethrow(exception);
        end
        refused=true;
    end
    result=test_true(result,['Full-observed-spin system workspace ' method{1}],refused,...
                     'Even one task must reserve full physical construction as well as bath propagation.');

    % A ceiling admitting one but not two workspaces must retain the high-budget result
    large_model.tols.cce_mem=resident+1.5*full_workspace;
    low=cce_schedule(large_model,large_clusters,large_ensemble,large_schedule,large);
    result=test_close(result,['Full-observed-spin system serial parity ' method{1}],...
                      low{1}(end).value,high{1}(end).value,0,0,...
                      'The memory ceiling reduces concurrency without changing the result.');

end

% Noise still constructs a full retained Liouville generator before bath restriction
noise_request=large;
noise_request.bath_jumps=struct('support',7,'operator',diag([1 -1])/2,'rate',1);
noise_model=large_model;
noise_model.tols.cce_mem=2^28;
noise_model.bas.formalism='zeeman-liouv';
refused=false;
try
    cce_schedule(noise_model,large_clusters,large_ensemble,large_schedule,...
                 noise_request);
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
large_model.bas.formalism='zeeman-liouv';
large_model.tols.cce_mem=2^28;
[subsystem,~]=cce_subsystem(large_model,liouville,7,struct(),[]);
[H,R,~]=cce_generators(subsystem,liouville,[]);
reference=-1i*(kron(eye(2),H{1})-kron(H{2}.',eye(2)))+R;
subsystem.tols.cce_mem=10000;
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
q=rmfield(p,{'bath_samples','bath_seed'});
observable_sys=spin_system;
observable_sys.tols.cce_mem=5000000;
q.rho0=eye(2)/2;
q.output='observable';
q.coil=repmat({eye(2)},1,200);
q.times=linspace(0,1,301);
q.bath_jumps=struct([]);
for method={'gcce','megcce'}
    q.cce.approximation='generalized';
    observable_sys.bas.formalism='zeeman-hilb';
    if strcmp(method{1},'megcce')
        observable_sys.bas.formalism='zeeman-liouv';
    end
    refused=false;
    try
        clexpn(observable_sys,@free_evolution,q,'labframe');
    catch exception
        refused=strcmp(exception.identifier,'cce_schedule:memory');
    end
    result=test_true(result,['observable column memory ' method{1}],refused,...
                     'The resident budget for four time-by-operator records exceeds five million bytes.');
    fprintf('R7 observable memory %s rejection %d raw minimum %d budget %d\n',...
            method{1},refused,4*25*numel(q.times)*numel(q.coil),observable_sys.tols.cce_mem);
end

% Retained ME manifold, spin, and jump operators accumulate across sampled tasks
for method={'mecce','megcce'}
    retained.observed=1;
    retained.times=[0 1e-3];
    retained.orientation=[0 0 0];
    retained.transition=struct('spins',1,'alpha',[1;0],'beta',[0;1]);
    retained.cce=struct('bath_spins',2,'order',1,'approximation','projected',...
                        'expansion','standard','interaction_cutoff',Inf,...
                        'meanfield','none','averaging','after_expansion',...
                        'zero_policy','mask');
    retained.bath_samples=40;
    retained.bath_seed=31;
    retained_sys=spin_system;
    retained_sys.tols.cce_atol=1e-12;
    retained_sys.tols.cce_rtol=1e-10;
    retained_sys.tols.cce_mem=1e8;
    retained.bath_jumps=repmat(struct('support',2,'operator',[0 1;0 0],'rate',2),1,8);
    if strcmp(method{1},'megcce')
        retained.cce.approximation='generalized';
        retained.output='density';
        retained.rho0=ones(2)/2;
    end
    retained=cce_options(retained_sys,retained,@free_evolution,'labframe');
    [~,retained_data]=clexpn(retained_sys,@free_evolution,retained,'labframe');
    retained_model_p=cce_options(retained_sys,retained,@cce_compile,'labframe');
    retained_model=cce_compile(retained_sys,retained_model_p,'labframe');
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
    retained_model.tols.cce_mem=resident+workspace+payload/2;
    refused=false;
    try
        cce_schedule(retained_model,retained_clusters,retained_states,...
                     retained_schedule,retained);
    catch exception
        refused=strcmp(exception.identifier,'cce_schedule:memory');
    end
    result=test_true(result,[method{1} ' retained metadata memory'],refused,...
                     'Every completed task retains its manifold, spin and jump operator arrays.');
    fprintf('%s retained operator payload %.0f bytes; cap %.0f; refused %d.\n',...
            method{1},payload,retained_model.tols.cce_mem,refused);
end

% Print the pilot timings without imposing an unjustified speedup target
fprintf('CCE pilot serial %.6g s; one worker %.6g s; two workers %.6g s; discrepancy %.3g.\n',...
        serial_time,one_time,two_time,max(abs(serial-two),[],'all'));
delete(gcp('nocreate'));
sys.parallel={'processes',1};
create(sys,inter);

end


