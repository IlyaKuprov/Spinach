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

end

