% Atomic checkpoint identity, interruption and corruption acceptance test.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_checkpoint()

% Declare exact restart and persistent-data integrity checks
result=new_test_result('kernel/cce_checkpoint','CCE atomic restart',...
                       'Partial ledger resume, identity rejection and write integrity.');

% Construct a small analytic Ising bath with deterministic sampled states
sys.isotopes={'13C','13C','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3); inter.coupling.matrix{1,2}=diag([0 0 17]);
inter.coupling.matrix{1,3}=diag([0 0 29]);
spin_system=create(sys,inter);
p=struct('central_spins',1,'bath_spins',[2 3],'method','cce',...
         'expansion','standard','quantity','coherence','times',[0 0.001 0.004],...
         'order',2,'central',struct('alpha',[1;0],'beta',[0;1]),...
         'assumptions','labframe','orientation',[0 0 0],'interaction_cutoff',Inf,...
         'graph',struct('mode','complete','score','coupling_norm'),...
         'meanfield','none','sequence',[],'averaging','after_expansion',...
         'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                             'samples',4,'seed',193),...
         'numerics',struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^26));

% Save a complete calculation and simulate a durable partial interruption
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
filename=fullfile(folder,'checkpoint.mat');
p.checkpoint=struct('file',filename,'resume',false);
p.parallel=struct('mode','serial','batch_size',3,'max_tasks',16);
uninterrupted=cce(spin_system,p); loaded=load(filename,'state'); state=loaded.state;
state.completed(5:end)=false; state.records(5:end)={[]};
cce_checkpoint(filename,'write',state);
p.checkpoint.resume=true; resumed=cce(spin_system,p);
result=test_true(result,'partial raw restart',isequaln(resumed.raw,uninterrupted.raw),...
                 'Only missing tasks are recomputed; all raw records are identical.');
result=test_close(result,'partial reduced restart',resumed.value,uninterrupted.value,0,0,...
                  'Restart repeats the same fixed-order reduction.');
loaded=load(filename,'state'); state=loaded.state;
result=test_true(result,'completed bitmap',all(state.completed),...
                 'All sixteen tasks are marked complete after resume.');

% Refuse cached records when any Spinach Hamiltonian tolerance changes
changed_system=spin_system; changed_system.tols.inter_cutoff=1e6; refused=false;
try
    cce(changed_system,p);
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
result=test_true(result,'Spinach tolerance checkpoint identity',refused,...
                 'A changed interaction cutoff must invalidate the physical checkpoint identity.');
q=rmfield(p,'checkpoint'); fresh=cce(changed_system,q);
tolerance_difference=max(abs(uninterrupted.value(:)-fresh.value(:)));
original_model=cce_prepare(spin_system,q); changed_model=cce_prepare(changed_system,q);
result=test_true(result,'Spinach tolerance result fingerprint',...
                 ~strcmp(original_model.source_hash,changed_model.source_hash),...
                 'The public model fingerprint must distinguish tolerance-dependent Hamiltonians.');
result=test_true(result,'Interaction cutoff changes propagation',tolerance_difference>1e-4,...
                 'The fixture must independently demonstrate that the changed tolerance affects results.');
fprintf('checkpoint tolerance observable difference %.12g\n',tolerance_difference);

% Refuse changed times, physical tensors, seed and cluster identity
changed=p; changed.times=[0 0.002 0.004]; refused=false;
try
    cce(spin_system,changed); catch exception, refused=strcmp(exception.identifier,'cce_checkpoint:identity'); end
result=test_true(result,'time identity',refused,'Time meshes are part of the restart identity.');
changed_model=spin_system; changed_model.inter.coupling.matrix{1,2}(3,3)=1.01*...
    changed_model.inter.coupling.matrix{1,2}(3,3); refused=false;
try
    cce(changed_model,p); catch exception, refused=strcmp(exception.identifier,'cce_checkpoint:identity'); end
result=test_true(result,'model identity',refused,'Physical tensors are hashed.');
changed=p; changed.bath_state.seed=194; refused=false;
try
    cce(spin_system,changed); catch exception, refused=strcmp(exception.identifier,'cce_checkpoint:identity'); end
result=test_true(result,'seed identity',refused,'Seed and explicit sampled states are hashed.');
changed=state; changed.identity.clusters='changed'; refused=false;
try
    cce_checkpoint(filename,'read',changed); catch exception, refused=strcmp(exception.identifier,'cce_checkpoint:identity'); end
result=test_true(result,'cluster identity',refused,'Cluster family and powers are hashed.');

% Check full-payload corruption and preserve good data on rejected writes
changed=state; changed.records{1}.value(1)=2; refused=false;
try
    cce_checkpoint(filename,'validate',changed); catch exception, refused=strcmp(exception.identifier,'cce_checkpoint:checksum'); end
result=test_true(result,'raw corruption',refused,'Checksum covers raw data, not just identity headers.');
changed=state; changed.records{1}.value(1)=2; refused=false;
try
    cce_checkpoint(filename,'write',changed);
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:completed');
end
result=test_true(result,'completed task immutability',refused,...
                 'Even matching identity cannot replace a different completed raw record.');
changed=state; changed.identity.model='different'; refused=false;
try
    cce_checkpoint(filename,'write',changed); catch exception, refused=strcmp(exception.identifier,'cce_checkpoint:identity'); end
loaded=load(filename,'state');
result=test_true(result,'failed replacement integrity',refused&&isequaln(loaded.state,state),...
                 'A failed replacement never modifies the known-good checkpoint.');
refused=false;
try
    cce_checkpoint(fullfile(folder,'missing','checkpoint.mat'),'write',state);
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:directory');
end
loaded=load(filename,'state');
result=test_true(result,'write failure integrity',refused&&isequaln(loaded.state,state),...
                 'A missing target directory cannot damage the existing checkpoint.');

% Force a real save failure with a read-only directory around a good file
[status,message]=fileattrib(folder,'-w','u'); %#ok<FILEATTRIB>
if ~status, error('test_cce_checkpoint:permissions','%s',message); end
refused=false;
try
    cce_checkpoint(filename,'write',state); catch
    refused=true;
end
[status,message]=fileattrib(folder,'+w','u'); %#ok<FILEATTRIB>
if ~status, error('test_cce_checkpoint:permissions','%s',message); end
loaded=load(filename,'state');
result=test_true(result,'actual save failure integrity',refused&&isequaln(loaded.state,state),...
                 'A real filesystem write failure leaves the known-good target byte content intact.');

% Verify execution settings do not falsely change physical restart identity
p.parallel.batch_size=1; p.parallel.mode='clusters'; changed=cce(spin_system,p);
result=test_close(result,'execution-neutral restart',changed.value,uninterrupted.value,0,0,...
                  'Worker batching and serial/parallel selection are not physics.');
fprintf('CCE checkpoint: %d tasks, 4 preserved, 12 resumed; discrepancy %.3g.\n',...
        numel(state.completed),max(abs(resumed.value-uninterrupted.value),[],'all'));

% Reject unsupported checkpoint ownership before direct evaluator work
plain=rmfield(p,{'checkpoint','parallel'}); plain.bath_state=struct('mode','mixed_product');
model=cce_prepare(spin_system,plain); ensemble=cce_states(model,plain);
schedule=cce_sequence(model,plain); cluster=struct('id',1,'sites',2);
for method={'cce','gcce','mecce','megcce','noise'}
    direct=plain; direct.method=method{1};
    if any(strcmp(direct.method,{'gcce','megcce'}))
        direct.central.rho=ones(2)/2; direct.normalized=false;
    end
    if any(strcmp(direct.method,{'mecce','megcce'}))
        direct.jumps=struct('support',2,'rate',2,'operator',[0 1;0 0]);
        evaluator=@cce_master;
    elseif strcmp(direct.method,'noise')
        direct.quantity='autocorrelation';
        direct.noise=struct('ordering','ordered','connected',false,...
                            'terms',struct('support',2,'operator',diag([1 -1])/2));
        evaluator=@cce_noise;
    else
        evaluator=@cce_cluster;
    end
    direct.checkpoint=struct('file',fullfile(folder,'direct.mat'),'resume',true);
    refused=false;
    try
        evaluator(model,cluster,ensemble.contexts(1),schedule,direct);
    catch exception
        refused=strcmp(exception.identifier,'cce:checkpoint');
    end
    result=test_true(result,['direct ' method{1} ' checkpoint rejection'],refused,...
                     'Only the scheduler owns checkpoint records; direct evaluators must not ignore them.');
end

% Reject interlaced checkpoint writes and resumes at both public entry points
interlaced=p; interlaced.bath_state.mode='interlaced';
interlaced.interlaced=struct('max_states',16,'reference_state',1);
model=cce_prepare(spin_system,rmfield(interlaced,'checkpoint'));
ensemble=cce_states(model,rmfield(interlaced,'checkpoint'));
schedule=cce_sequence(model,rmfield(interlaced,'checkpoint'));
for resume=[false true]
    interlaced.checkpoint=struct('file',fullfile(folder,'interlaced.mat'),'resume',resume);
    for direct=[false true]
        refused=false;
        try
            if direct
                cce_interlace(model,cluster,ensemble,schedule,interlaced);
            else
                cce(spin_system,interlaced);
            end
        catch exception
            refused=strcmp(exception.identifier,'cce:checkpoint');
        end
        result=test_true(result,sprintf('interlaced checkpoint resume=%d direct=%d',resume,direct),...
                         refused&&~isfile(interlaced.checkpoint.file),...
                         'Interlaced runs have no checkpoint ledger and must reject write and resume requests.');
    end
end

% Reject malformed checkpoint requests before any physical work
for specification={struct('file',''),struct('file','','resume',true),struct('resume',false)}
    request=plain; request.checkpoint=specification{1}; refused=false;
    try
        cce(spin_system,request);
    catch exception
        refused=strcmp(exception.identifier,'cce:checkpoint');
    end
    result=test_true(result,'checkpoint filename required',refused,...
                     'An empty or absent filename must not silently disable a requested checkpoint.');
end

% Keep adaptive convergence explicitly incompatible with a shared checkpoint
adaptive=p; adaptive.convergence=struct('ladders',struct('field','order','values',{{1 2}}),...
                                       'atol',1e-8,'rtol',1e-6,'max_work',3);
refused=false;
try
    cce_converge(spin_system,adaptive);
catch exception
    refused=strcmp(exception.identifier,'cce_converge:checkpoint');
end
result=test_true(result,'adaptive checkpoint rejection',refused,...
                 'Independent convergence configurations cannot share a single checkpoint ledger.');

% Exercise restart identity across all scheduler evaluator and control paths
for method={'cce','gcce','mecce','megcce','noise','controls'}
    request=plain;
    if strcmp(method{1},'controls')
        request.method='gcce'; request.central.rho=ones(2)/2;
        request.sequence=struct('type','finite','target',1,'placement','absolute',...
                                'at',0,'duration',0.0005,'axis','x','amplitude',7,'detuning',3);
        request.times=[0 0.002 0.004]; request.normalized=false;
    else
        request.method=method{1};
    end
    if any(strcmp(request.method,{'gcce','megcce'}))
        request.central.rho=ones(2)/2; request.normalized=false;
    end
    if any(strcmp(request.method,{'mecce','megcce'}))
        request.jumps=struct('support',2,'rate',2,'operator',[0 1;0 0]);
    elseif strcmp(request.method,'noise')
        request.quantity='autocorrelation';
        request.noise=struct('ordering','ordered','connected',false,...
                             'terms',struct('support',2,'operator',diag([1 -1])/2));
    end
    request.checkpoint=struct('file',fullfile(folder,[method{1} '.mat']),'resume',false);
    request.parallel=struct('mode','serial','batch_size',2,'max_tasks',4);
    complete=cce(spin_system,request); saved=load(request.checkpoint.file,'state');
    partial=saved.state; partial.completed(2:end)=false; partial.records(2:end)={[]};
    cce_checkpoint(request.checkpoint.file,'write',partial);
    request.checkpoint.resume=true; request.parallel.mode='clusters';
    restarted=cce(spin_system,request);
    result=test_close(result,[method{1} ' partial serial parallel restart'],...
                      restarted.value,complete.value,0,0,...
                      'Changing worker scheduling must preserve the fixed-order resumed observable exactly.');
    changed=request; changed.times=2*request.times; refused=false;
    try
        cce(spin_system,changed);
    catch exception
        refused=strcmp(exception.identifier,'cce_checkpoint:identity');
    end
    result=test_true(result,[method{1} ' changed configuration rejection'],refused,...
                     'A checkpoint from another time grid cannot supply cached physical results.');
    % Reject changes to evaluator-specific physics as well as the time grid
    changed=request;
    if strcmp(method{1},'controls')
        changed.sequence.amplitude=2*request.sequence.amplitude;
    elseif any(strcmp(method{1},{'mecce','megcce'}))
        changed.jumps.rate=2*request.jumps.rate;
    elseif strcmp(method{1},'noise')
        changed.noise.terms.operator=2*request.noise.terms.operator;
    else
        changed.orientation=[0.2 0.3 0.4];
    end
    refused=false;
    try
        cce(spin_system,changed);
    catch exception
        refused=strcmp(exception.identifier,'cce_checkpoint:identity');
    end
    result=test_true(result,[method{1} ' changed physics rejection'],refused,...
                     'Each evaluator-specific physical option participates in checkpoint identity.');

    fprintf('R14 restart %s discrepancy %.17g raw identical %d\n',method{1},...
            max(abs(restarted.value-complete.value),[],'all'),isequaln(restarted.raw,complete.raw));
end

% Compare successive evaluations with cold references after physical changes
for method={'cce','mecce','noise','controls'}
    request=plain;
    if strcmp(method{1},'mecce')
        request.method='mecce';
        request.jumps=struct('support',2,'rate',2,'operator',[0 1;0 0]);
    elseif strcmp(method{1},'noise')
        request.method='noise'; request.quantity='autocorrelation';
        request.noise=struct('ordering','ordered','connected',false,...
                             'terms',struct('support',2,'operator',diag([1 -1])/2));
    elseif strcmp(method{1},'controls')
        request.method='gcce'; request.central.rho=ones(2)/2; request.normalized=false;
        request.sequence=struct('type','finite','target',1,'placement','absolute',...
                                'at',0,'duration',0.0005,'axis','x','amplitude',7,'detuning',3);
        request.times=[0 0.002 0.004];
    end
    original=cce(spin_system,request); changed_system=spin_system;
    changed_system.inter.coupling.matrix{1,2}=2*spin_system.inter.coupling.matrix{1,2};
    warm=cce(changed_system,request); clear cce_hamiltonian;
    cold=cce(changed_system,request);
    result=test_close(result,[method{1} ' successive cleared-function identity'],warm.value,cold.value,0,0,...
                      'Changing the physical model must not reuse any previous physical result.');
    result=test_true(result,[method{1} ' physical model fingerprint'],...
                     ~strcmp(original.model_hash,warm.model_hash),...
                     'Different physical interactions require distinct model fingerprints.');
end

end


