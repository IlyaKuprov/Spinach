% Test: kernel/cce_checkpoint | CCE atomic restart
% Atomic checkpoint identity, interruption, and corruption acceptance test. Syntax:
%
%                    result=test_cce_checkpoint()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_checkpoint()

% Declare exact restart and persistent-data integrity checks
result=new_test_result('kernel/cce_checkpoint','CCE atomic restart',...
                       'Partial ledger resume, identity rejection and write integrity.');

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

% Save a complete calculation and simulate a durable partial interruption
folder=tempname;
mkdir(folder);
cleanup=onCleanup(@() rmdir(folder,'s'));
filename=fullfile(folder,'checkpoint.mat');
p.cce.checkpoint=struct('file',filename,'resume',false);
p.cce.parallel=struct('mode','serial','batch_size',3,'max_tasks',16);
[uninterrupted,uninterrupted_data]=clexpn(spin_system,@free_evolution,p,'labframe');
loaded=load(filename,'state');
state=loaded.state;
state.completed(5:end)=false;
state.records(5:end)={[]};
cce_checkpoint(filename,'write',state);
p.cce.checkpoint.resume=true;
[resumed,resumed_data]=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_true(result,'partial raw restart',isequaln(resumed_data.raw,uninterrupted_data.raw),...
                 'Only missing tasks are recomputed; all raw records are identical.');
result=test_close(result,'partial reduced restart',resumed,uninterrupted,0,0,...
                  'Restart repeats the same fixed-order reduction.');
loaded=load(filename,'state');
state=loaded.state;
result=test_true(result,'completed bitmap',all(state.completed),...
                 'All sixteen tasks are marked complete after resume.');

% Only a retained physical orientation contributes to its grid checkpoint
quadrature=rmfield(p,'orientation');
quadrature.grid=struct('alphas',[1 0],'betas',[1 0],'gammas',[1 0],'weights',[0 3]);
[written,grid_data]=clexpn(spin_system,@free_evolution,quadrature,'labframe');
result=test_close(result,'Single-crystal to grid restart',written,uninterrupted,0,0,...
                  'A single retained Euler point can resume its direct crystal checkpoint.');
quadrature.grid.betas(1)=0.2;
quadrature.grid.weights(2)=7;
replayed=clexpn(spin_system,@free_evolution,quadrature,'labframe');
result=test_true(result,'Single-node checkpoint records',isequal(grid_data.grid_indices,2)&&...
                 isequaln(grid_data.orientations{1}.raw,uninterrupted_data.raw),...
                 'The grid node and direct crystal have identical physical records.');
result=test_close(result,'Single-node checkpoint replay',replayed,written,0,0,...
                  'Changing a discarded node and the relative weight preserves the grid ledger.');
quadrature.grid.betas(2)=0.2;
refused=false;
try
    clexpn(spin_system,@free_evolution,quadrature,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
result=test_true(result,'Retained orientation checkpoint identity',refused,...
                 'A different physical rotation must not resume the previous orientation.');

% Resource availability must not invalidate an otherwise identical physical ledger
resized=p;
resized.cce.numerics.max_memory_bytes=2*p.cce.numerics.max_memory_bytes;
[resource_resumed,resource_data]=clexpn(spin_system,@free_evolution,resized,'labframe');
result=test_true(result,'Resource-independent checkpoint identity',...
                 isequaln(resource_data.raw,uninterrupted_data.raw),...
                 'Changing only the memory allowance preserves cached physical records.');
result=test_close(result,'Resource-independent checkpoint replay',resource_resumed,uninterrupted,0,0,...
                  'The new resource limit is checked before the same physical ledger is resumed.');

% Keep adaptive pilot and held-out evaluations separate from checkpoint ownership
adaptive=p;
adaptive.cce=rmfield(adaptive.cce,'checkpoint');
adaptive.cce.order=1;
adaptive.cce.expansion='partition';
pilot=p;
pilot.cce=rmfield(pilot.cce,'checkpoint');
pilot.cce.order=0;
[~,pilot_data]=clexpn(spin_system,@free_evolution,pilot,'labframe');
adaptive_states=pilot_data.ensemble;
adaptive.cce.partition=struct('mode','adaptive','max_dimension',4,...
                              'patch_atol',1e-10,...
                              'patch_rtol',1e-10,'max_rounds',2,'max_probes',100,'patches',{{[2 3]}},...
                              'pilot_configs',adaptive_states.spin_config,'heldout_configs',adaptive_states.spin_config);
reference=clexpn(spin_system,@free_evolution,adaptive,'labframe');
accepted=false;
discrepancy=Inf;
adaptive.cce.checkpoint=struct('file',fullfile(folder,'adaptive-partition.mat'),'resume',false);
try
    [fresh,fresh_data]=clexpn(spin_system,@free_evolution,adaptive,'labframe');
    adaptive.cce.checkpoint.resume=true;
    [replayed,replayed_data]=clexpn(spin_system,@free_evolution,adaptive,'labframe');
    accepted=isequaln(fresh_data.raw,replayed_data.raw);
    discrepancy=max(abs([fresh-reference;replayed-reference]),[],'all');
catch exception
    if ~strcmp(exception.identifier,'cce_options:checkpoint')
        rethrow(exception);
    end
end
result=test_true(result,'Adaptive partition checkpoint ownership',accepted,...
                 'Pilot and held-out patches must not inherit the production checkpoint.');
result=test_close(result,'Adaptive partition checkpoint replay',discrepancy,0,1e-13,0,...
                  'Fresh and resumed adaptive partitions must agree with the checkpoint-free result.');
fprintf('Adaptive checkpoint raw equality %d; discrepancy %.3g.\n',accepted,discrepancy);

% Refuse cached records when any Spinach Hamiltonian tolerance changes
changed_system=spin_system;
changed_system.tols.inter_cutoff=1e6;
refused=false;
try
    clexpn(changed_system,@free_evolution,p,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
result=test_true(result,'Spinach tolerance checkpoint identity',refused,...
                 'A changed interaction cutoff must invalidate the physical checkpoint identity.');
q=p;
q.cce=rmfield(q.cce,'checkpoint');
[fresh,~]=clexpn(changed_system,@free_evolution,q,'labframe');
tolerance_difference=max(abs(uninterrupted(:)-fresh(:)));
original_model=cce_compile(spin_system,q,'labframe');
changed_model=cce_compile(changed_system,q,'labframe');
result=test_true(result,'Spinach tolerance result fingerprint',...
                 ~strcmp(original_model.source_hash,changed_model.source_hash),...
                 'The public model fingerprint must distinguish tolerance-dependent Hamiltonians.');
result=test_true(result,'Interaction cutoff changes propagation',tolerance_difference>1e-4,...
                 'The fixture must independently demonstrate that the changed tolerance affects results.');
fprintf('checkpoint tolerance observable difference %.12g\n',tolerance_difference);

% Refuse changed times, physical tensors, seed, and cluster identity
changed=p;
changed.times=[0 0.002 0.004];
refused=false;
try
    clexpn(spin_system,@free_evolution,changed,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
result=test_true(result,'time identity',refused,'Time meshes are part of the restart identity.');
changed_model=spin_system;
changed_model.inter.coupling.matrix{1,2}(3,3)=1.01*...
    changed_model.inter.coupling.matrix{1,2}(3,3);
refused=false;
try
    clexpn(changed_model,@free_evolution,p,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
result=test_true(result,'model identity',refused,'Physical tensors are hashed.');
changed=p;
changed.cce.bath_state.seed=194;
refused=false;
try
    clexpn(spin_system,@free_evolution,changed,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
result=test_true(result,'seed identity',refused,'Seed and explicit sampled states are hashed.');
changed=state;
changed.identity.clusters='changed';
refused=false;
try
    cce_checkpoint(filename,'read',changed);
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
result=test_true(result,'cluster identity',refused,'Cluster family and powers are hashed.');

% Check full-payload corruption and preserve good data on rejected writes
changed=state;
changed.records{1}.value(1)=2;
refused=false;
try
    cce_checkpoint(filename,'validate',changed);
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:checksum');
end
result=test_true(result,'raw corruption',refused,'Checksum covers raw data, not just identity headers.');
changed=state;
changed.records{1}.value(1)=2;
refused=false;
try
    cce_checkpoint(filename,'write',changed);
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:completed');
end
result=test_true(result,'completed task immutability',refused,...
                 'Even matching identity cannot replace a different completed raw record.');
changed=state;
changed.identity.model='different';
refused=false;
try
    cce_checkpoint(filename,'write',changed);
catch exception
    refused=strcmp(exception.identifier,'cce_checkpoint:identity');
end
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

% Deny new-file creation around a known-good checkpoint
if ispc
    [status,message]=system(sprintf('icacls "%s" /deny "*S-1-1-0:(WD,AD)"',folder));
    status=(status==0);
else
    [status,message]=fileattrib(folder,'-w','u'); %#ok<FILEATTRIB>
end
if ~status
    error('test_cce_checkpoint:permissions','%s',message);
end
refused=false;
try
    cce_checkpoint(filename,'write',state);
catch
    refused=true;
end
if ispc
    [status,message]=system(sprintf('icacls "%s" /remove:d "*S-1-1-0"',folder));
    status=(status==0);
else
    [status,message]=fileattrib(folder,'+w','u'); %#ok<FILEATTRIB>
end
if ~status
    error('test_cce_checkpoint:permissions','%s',message);
end
loaded=load(filename,'state');
result=test_true(result,'actual save failure integrity',refused&&isequaln(loaded.state,state),...
                 'A real filesystem write failure leaves the known-good target byte content intact.');

% Verify execution settings do not falsely change physical restart identity
p.cce.parallel.batch_size=1;
p.cce.parallel.mode='clusters';
changed=clexpn(spin_system,@free_evolution,p,'labframe');
result=test_close(result,'execution-neutral restart',changed,uninterrupted,0,0,...
                  'Worker batching and serial/parallel selection are not physics.');
fprintf('CCE checkpoint: %d tasks, 4 preserved, 12 resumed; discrepancy %.3g.\n',...
        numel(state.completed),max(abs(resumed-uninterrupted),[],'all'));

% Reject unsupported checkpoint ownership before direct evaluator work
plain=p;
plain.cce=rmfield(plain.cce,{'checkpoint','parallel'});
plain.cce.bath_state=struct('mode','mixed_product');
model=cce_compile(spin_system,plain,'labframe');
ensemble=uninterrupted_data.ensemble;
cluster=struct('id',1,'sites',2);
for method={'cce','gcce','mecce','megcce','noise'}
    direct=plain;
    pulse_sequence=@free_evolution;
    if any(strcmp(method{1},{'gcce','megcce'}))
        direct.cce.approximation='generalized';
        direct.rho0=ones(2)/2;
        direct.normalized=false;
    end
    if any(strcmp(method{1},{'mecce','megcce'}))
        direct.bath_jumps=struct('support',2,'rate',2,'operator',[0 1;0 0]);
        direct.cce.basis.formalism='zeeman-liouv';
    elseif strcmp(method{1},'noise')
        pulse_sequence=@noise_correlation;
        direct.cce.approximation='generalized';
        direct.noise=struct('ordering','ordered','connected',false,'manifold','bath_only',...
                            'descriptors',struct('spins',2,'operators','Lz','coefficient',1));
    end
    schedule=struct('pulse_sequence',pulse_sequence,'contract',cce_capabilities(pulse_sequence));
    evaluator=schedule.contract.evaluator;
    direct.cce.checkpoint=struct('file',fullfile(folder,'direct.mat'),'resume',true);
    refused=false;
    try
        evaluator(model,cluster,ensemble.spin_config(1),schedule,direct);
    catch exception
        refused=strcmp(exception.identifier,'cce_options:checkpoint');
    end
    result=test_true(result,['direct ' method{1} ' checkpoint rejection'],refused,...
                     'Only the scheduler owns checkpoint records; direct evaluators must not ignore them.');
end

% Reject interlaced checkpoint writes and resumes at both public entry points
interlaced=p;
interlaced.cce.bath_state.mode='interlaced';
interlaced.cce.interlaced=struct('max_states',16,'reference_state',1);
model=cce_compile(spin_system,plain,'labframe');
ensemble=uninterrupted_data.ensemble;
schedule=struct('pulse_sequence',@free_evolution,'contract',cce_capabilities(@free_evolution));
for resume=[false true]
    interlaced.cce.checkpoint=struct('file',fullfile(folder,'interlaced.mat'),'resume',resume);
    for direct=[false true]
        refused=false;
        try
            if direct
                cce_interlace(model,cluster,ensemble,schedule,interlaced,schedule.contract.evaluator);
            else
                clexpn(spin_system,@free_evolution,interlaced,'labframe');
            end
        catch exception
            refused=strcmp(exception.identifier,'cce_options:checkpoint');
        end
        result=test_true(result,sprintf('interlaced checkpoint resume=%d direct=%d',resume,direct),...
                         refused&&~isfile(interlaced.cce.checkpoint.file),...
                         'Interlaced runs have no checkpoint ledger and must reject write and resume requests.');
    end
end

% Reject malformed checkpoint requests before any physical work
for specification={struct('file',''),struct('file','','resume',true),struct('resume',false)}
    request=plain;
    request.cce.checkpoint=specification{1};
    refused=false;
    try
        clexpn(spin_system,@free_evolution,request,'labframe');
    catch exception
        refused=strcmp(exception.identifier,'cce_options:checkpoint');
    end
    result=test_true(result,'checkpoint filename required',refused,...
                     'An empty or absent filename must not silently disable a requested checkpoint.');
end

% Keep adaptive convergence explicitly incompatible with a shared checkpoint
adaptive=p;
adaptive.cce.convergence=struct('ladders',struct('field','order','values',{{1 2}}),...
                                'atol',1e-8,'rtol',1e-6,'max_work',3);
refused=false;
try
    clexpn(spin_system,@free_evolution,adaptive,'labframe');
catch exception
    refused=strcmp(exception.identifier,'cce_options:checkpoint');
end
result=test_true(result,'adaptive checkpoint rejection',refused,...
                 'Independent convergence configurations cannot share a single checkpoint ledger.');

% Exercise restart identity across all scheduler evaluator and control paths
for method={'cce','gcce','mecce','megcce','noise','controls'}
    request=plain;
    pulse_sequence=@free_evolution;
    if any(strcmp(method{1},{'gcce','megcce','controls'}))
        request.cce.approximation='generalized';
        request.rho0=ones(2)/2;
        request.normalized=false;
    end
    if strcmp(method{1},'controls')
        pulse_sequence=@finite_drive;
        request.drive=struct('spins',1,'duration',0.0005,...
                             'axis','x','amplitude',7,'detuning',3);
        request.times=[0 0.002 0.004];
    elseif any(strcmp(method{1},{'mecce','megcce'}))
        request.bath_jumps=struct('support',2,'rate',2,'operator',[0 1;0 0]);
        request.cce.basis.formalism='zeeman-liouv';
    elseif strcmp(method{1},'noise')
        pulse_sequence=@noise_correlation;
        request.cce.approximation='generalized';
        request.noise=struct('ordering','ordered','connected',false,'manifold','bath_only',...
                             'descriptors',struct('spins',2,'operators','Lz','coefficient',1));
    end
    request.cce.checkpoint=struct('file',fullfile(folder,[method{1} '.mat']),'resume',false);
    request.cce.parallel=struct('mode','serial','batch_size',2,'max_tasks',4);
    [complete,complete_data]=clexpn(spin_system,pulse_sequence,request,'labframe');
    saved=load(request.cce.checkpoint.file,'state');
    partial=saved.state;
    partial.completed(2:end)=false;
    partial.records(2:end)={[]};
    cce_checkpoint(request.cce.checkpoint.file,'write',partial);
    request.cce.checkpoint.resume=true;
    request.cce.parallel.mode='clusters';
    [restarted,restarted_data]=clexpn(spin_system,pulse_sequence,request,'labframe');
    result=test_close(result,[method{1} ' partial serial parallel restart'],...
                      restarted,complete,0,0,...
                      'Changing worker scheduling must preserve the fixed-order resumed observable exactly.');
    result=test_true(result,[method{1} ' resumed raw identity'],isequaln(restarted_data.raw,complete_data.raw),...
                     'Resumed and uninterrupted raw records must be identical on every evaluator path.');
    changed=request;
    changed.times=2*request.times;
    refused=false;
    try
        clexpn(spin_system,pulse_sequence,changed,'labframe');
    catch exception
        refused=strcmp(exception.identifier,'cce_checkpoint:identity');
    end
    result=test_true(result,[method{1} ' changed configuration rejection'],refused,...
                     'A checkpoint from another time grid cannot supply cached physical results.');

    % Reject changes to evaluator-specific physics as well as the time grid
    changed=request;
    if strcmp(method{1},'controls')
        changed.drive.amplitude=2*request.drive.amplitude;
    elseif any(strcmp(method{1},{'mecce','megcce'}))
        changed.bath_jumps.rate=2*request.bath_jumps.rate;
    elseif strcmp(method{1},'noise')
        changed.noise.descriptors.coefficient=2*request.noise.descriptors.coefficient;
    else
        changed.orientation=[0.2 0.3 0.4];
    end
    refused=false;
    try
        clexpn(spin_system,pulse_sequence,changed,'labframe');
    catch exception
        refused=strcmp(exception.identifier,'cce_checkpoint:identity');
    end
    result=test_true(result,[method{1} ' changed physics rejection'],refused,...
                     'Each evaluator-specific physical option participates in checkpoint identity.');

    fprintf('R14 restart %s discrepancy %.17g raw identical %d\n',method{1},...
            max(abs(restarted-complete),[],'all'),isequaln(restarted_data.raw,complete_data.raw));
end

% Compare successive evaluations with cold references after physical changes
for method={'cce','mecce','noise','controls'}
    request=plain;
    pulse_sequence=@free_evolution;
    if any(strcmp(method{1},{'gcce','megcce','controls'}))
        request.cce.approximation='generalized';
        request.rho0=ones(2)/2;
        request.normalized=false;
    end
    if strcmp(method{1},'controls')
        pulse_sequence=@finite_drive;
        request.drive=struct('spins',1,'duration',0.0005,...
                             'axis','x','amplitude',7,'detuning',3);
        request.times=[0 0.002 0.004];
    elseif any(strcmp(method{1},{'mecce','megcce'}))
        request.bath_jumps=struct('support',2,'rate',2,'operator',[0 1;0 0]);
        request.cce.basis.formalism='zeeman-liouv';
    elseif strcmp(method{1},'noise')
        pulse_sequence=@noise_correlation;
        request.cce.approximation='generalized';
        request.noise=struct('ordering','ordered','connected',false,'manifold','bath_only',...
                             'descriptors',struct('spins',2,'operators','Lz','coefficient',1));
    end
    original=clexpn(spin_system,pulse_sequence,request,'labframe'); %#ok<NASGU>
    changed_system=spin_system;
    changed_system.inter.coupling.matrix{1,2}=2*spin_system.inter.coupling.matrix{1,2};
    warm=clexpn(changed_system,pulse_sequence,request,'labframe');
    clear cce_hamiltonian;
    cold=clexpn(changed_system,pulse_sequence,request,'labframe');
    result=test_close(result,[method{1} ' successive cleared-function identity'],warm,cold,0,0,...
                      'Changing the physical model must not reuse any previous physical result.');
    before_model=cce_compile(spin_system,request,'labframe');
    after_model=cce_compile(changed_system,request,'labframe');
    result=test_true(result,[method{1} ' physical model fingerprint'],...
                     ~strcmp(before_model.source_hash,after_model.source_hash),...
                     'Different physical interactions require distinct model fingerprints.');
end

end


