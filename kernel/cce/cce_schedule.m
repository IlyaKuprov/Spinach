% Evaluate CCE state/cluster tasks in deterministic indexed batches. Syntax:
%
%     raw=cce_schedule(spin_system,clusters,ensemble,pulse_sequence,parameters)
%
% raw is a cell per input spin_config, containing the empty-cluster record
% followed by cluster records in input order. Tasks execute in ascending
% configuration/state/cluster ID order; no reduction occurs on workers.
% Tasks run on the pool that create started from sys.parallel, and in the
% client when that pool has one worker or the caller is itself a worker.
% sys.tols.cce_mem bounds resident raw data plus simultaneous per-worker
% dense storage; the worker count is reduced until that bound is met. The
% memory preflight includes retained master generators and full physical
% construction before projection. parameters.cce.checkpoint.file/resume
% enable atomic restart, written after every batch of one task per worker.
%
% Parameters:
%
%   spin_system, clusters, ensemble, pulse_sequence, parameters - prepared
%   objects and the experiment function handle. Optional
%   parameters.cce.checkpoint.file/resume.
%
% Outputs:
%
%   raw - cell per input spin_config, each containing the baseline record then
%   cluster records in input order.
%
% See also cce_checkpoint, cce_evaluate, clexpn.
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_schedule(spin_system,clusters,ensemble,pulse_sequence,parameters)

% Check consistency
grumble(spin_system,clusters,ensemble,pulse_sequence,parameters);

% Build deterministic tasks from resolved native controls
options=parameters.cce;
ns=numel(ensemble.spin_config);
nc=numel(clusters)+1;
components=1;
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    components=prod(spin_system.comp.mults(parameters.observed))^2;
elseif isfield(parameters,'coil')&&iscell(parameters.coil)
    components=numel(parameters.coil);
end
nt=numel(parameters.times);
dc=prod(spin_system.comp.mults(parameters.observed));
metadata=whos('spin_system','clusters','ensemble','parameters');
resident=3*ns*nc*(nt*(25*components+16*dc^2)+1024)+sum([metadata.bytes]);
noise=isequal(pulse_sequence,@noise_correlation);
master=strcmp(spin_system.bas.formalism,'zeeman-liouv')||...
       (isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps));
if noise
    resident=resident+3*ns*nc*(32*nt+16);
end

% Count retained experiment-owned control intervals before dispatch
segment_count=0;
if isequal(pulse_sequence,@finite_drive)||isequal(pulse_sequence,@refocusing_train)
    events=cce_pulse_events(parameters);
    for t=parameters.times
        starts=[events.at].';
        stops=starts;
        for j=1:numel(events)
            if strcmp(events(j).placement,'fraction')
                starts(j)=starts(j)*t;
            end
            stops(j)=starts(j);
            if strcmp(events(j).type,'finite')
                stops(j)=starts(j)+events(j).duration;
            end
        end
        segment_count=segment_count+numel(unique([0;t;starts(starts>0&starts<t);stops(stops>0&stops<t)]));
    end
end

% Reserve full physical construction and every simultaneous propagation workspace
peak=0;
for k=0:numel(clusters)
    sites=[];
    if k>0
        sites=clusters(k).sites;
    end
    ids=[parameters.observed sites];
    db=prod(spin_system.comp.mults(sites));
    dimension=dc*db;
    copies=12+3*numel(ids);
    if ~spin_system.cce.direct_supported||strcmp(parameters.cce.backend,'stock')
        ranks=max([2 cellfun(@numel,spin_system.inter.giant.coeff).']);
        copies=copies+sum((2*(1:ranks)+1).^2);
    end
    workspace=16*copies*dimension^2;
    propagated=dimension;
    if (strcmp(options.approximation,'projected')&&~noise)||(noise&&~master)
        propagated=db;
    end
    if master

        % Count manifold, spin, and included jump arrays retained by every task
        jumps=0;
        labels=0;
        if isfield(parameters,'bath_jumps')
            for jump=parameters.bath_jumps
                if all(ismember(jump.support,ids))
                    jumps=jumps+1;
                    labels=labels+numel(jump.support)+1;
                end
            end
        end
        retained=16*(3*numel(ids)*dimension^2+2*(1+jumps)*propagated^2)+...
                 8*(2*numel(ids)+labels);
        resident=resident+3*ns*retained;
        workspace=workspace+retained;
        propagated=propagated^2;

    end
    workspace=workspace+16*24*propagated^2+segment_count*(48*dimension^2+24*numel(parameters.cce.bath_spins)+8);
    peak=max(peak,workspace);
end

% Use as many pool workers as the memory budget admits; none inside a worker
workers=min([max(1,poolsize()) ns*nc floor((spin_system.tols.cce_mem-resident)/peak)]);
if workers<1
    error('cce_schedule:memory','Raw ledger and one dense task exceed sys.tols.cce_mem.');
end

% Keep the full family out of conventional worker payloads
family=[];
if isfield(options,'interlaced')
    family=clusters;
end
task_clusters=[struct('id',uint64(0),'sites',[]);...
               struct('id',{clusters.id},'sites',{clusters.sites}).'];

% Hash physical and experiment inputs, not checkpoint filenames or resume switches
payload=struct('model',spin_system,'clusters',clusters,'ensemble',ensemble,...
               'parameters',parameters,'experiment',func2str(pulse_sequence));
for field={'checkpoint','retain_raw'}
    if isfield(payload.parameters.cce,field{1})
        payload.parameters.cce=rmfield(payload.parameters.cce,field{1});
    end
end

% The spin system is hashed from its full physical content, not its derived cache key
payload.model.cce=rmfield(payload.model.cce,'source_hash');
ids=zeros(ns*nc,3);
indices=zeros(ns*nc,2);
position=0;
for n=1:ns
    for k=1:nc
        position=position+1;
        id=0;
        if k>1
            id=clusters(k-1).id;
        end
        ids(position,:)=[ensemble.spin_config(n).configuration_id ensemble.spin_config(n).state_id id];
        indices(position,:)=[n k];
    end
end
[ids,order]=sortrows(ids,[1 2 3]);
indices=indices(order,:);
positions=zeros(size(order));
positions(order)=1:numel(order);
identity=struct();
if isfield(options,'checkpoint')

    % Available memory is a dispatch limit, not part of the physical calculation
    physical=payload;
    physical.model.tols=rmfield(physical.model.tols,'cce_mem');
    identity=struct('model',md5_hash(physical.model),'parameters',md5_hash(physical.parameters),...
                    'clusters',md5_hash(clusters),'schedule',md5_hash(payload.experiment),...
                    'samples',md5_hash(ensemble));

end
state=struct('version',1,'identity',identity,'payload',payload,...
             'task_ids',ids,'completed',false(ns*nc,1),'records',{cell(ns*nc,1)},...
             'checksum','');
filename='';
if isfield(options,'checkpoint')
    filename=options.checkpoint.file;
    if isfield(options.checkpoint,'resume')&&options.checkpoint.resume
        state.physical=physical;
        state=cce_checkpoint(filename,'read',state);
    elseif isfile(filename)
        error('cce_schedule:exists','Explicitly resume an existing checkpoint or choose another filename.');
    end
end

% Keep checkpoint ownership in the scheduler, not its experiment callback
if isfield(parameters.cce,'checkpoint')
    parameters.cce=rmfield(parameters.cce,'checkpoint');
end
pending=find(~state.completed);

% Write a checkpoint after every batch of one task per worker
batch_size=max(1,numel(pending));
if ~isempty(filename)
    batch_size=workers;
end
for b=1:batch_size:numel(pending)
    selected=pending(b:min(b+batch_size-1,numel(pending)));
    records=cell(numel(selected),1);
    tasks=indices(selected,:);
    local_clusters=task_clusters(tasks(:,2));
    configuration_ids=tasks(:,1);
    if workers<2
        for j=1:numel(selected)
            records{j}=native_task(spin_system,local_clusters(j),family,ensemble,pulse_sequence,parameters,configuration_ids(j));
        end
    else
        parfor (j=1:numel(selected),workers)
            records{j}=native_task(spin_system,local_clusters(j),family,ensemble,pulse_sequence,parameters,configuration_ids(j));
        end
    end
    state.records(selected)=records;
    state.completed(selected)=true;
    if ~isempty(filename)
        state=cce_checkpoint(filename,'write',state);
    end
end
raw=cell(1,ns);
for n=1:ns
    records=cell(1,nc);
    for k=1:nc
        position=positions((n-1)*nc+k);
        records{k}=state.records{position};
    end
    if isfield(options,'interlaced')
        records=records(~cellfun(@(item)strcmp(item.kind,'not_evaluated'),records));
    end
    raw{n}=[records{:}];
end

end

% Evaluate one cluster on the ordinary or the conditioned path
function raw=native_task(spin_system,cluster,family,ensemble,pulse_sequence,parameters,config_id)

% Keep conditioned tasks on the same physical family and configuration
if isfield(parameters.cce,'interlaced')
    parameters.cce.interlaced.reference_state=config_id;
    complete=any(arrayfun(@(item)isequal(sort(item.sites),sort(parameters.cce.bath_spins)),family));
    if complete&&~isequal(sort(cluster.sites),sort(parameters.cce.bath_spins))
        channels=1;
        if isfield(parameters,'coil')&&iscell(parameters.coil)
            channels=numel(parameters.coil);
        end
        raw=struct('cluster_id',cluster.id,'value',NaN(numel(parameters.times),channels),...
                   'abs_error_bound',Inf(numel(parameters.times),channels),...
                   'valid',false(numel(parameters.times),channels),...
                   'weight',ensemble.spin_config(config_id).weight,...
                   'condition_sites',{{}},'assignment_count',0,'reference',[],...
                   'reference_error',[],'full_bath',false,'kind','not_evaluated');
        if strcmp(spin_system.bas.formalism,'zeeman-liouv')||...
           (isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps))
            raw.generator=[];
        end
    else
        raw=cce_interlace(spin_system,cluster,family,ensemble,pulse_sequence,parameters);
    end
else
    raw=cce_evaluate(spin_system,cluster,ensemble.spin_config(config_id),pulse_sequence,parameters);
end
end

% Consistency enforcement
function grumble(spin_system,clusters,ensemble,pulse_sequence,parameters)
if ~isstruct(spin_system)||~isscalar(spin_system)||~isfield(spin_system,'cce')||...
   ~isstruct(clusters)||(~isempty(clusters)&&~all(isfield(clusters,{'id','sites'})))||...
   ~isstruct(ensemble)||~isfield(ensemble,'spin_config')||isempty(ensemble.spin_config)||...
   ~isa(pulse_sequence,'function_handle')||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'cce')
    error('cce_schedule:input','Native spin_system, ensemble, experiment handle and parameters are required.');
end
if isfield(parameters.cce,'interlaced')&&~isfield(parameters,'bath_samples')
    error('cce_schedule:interlaced','Interlaced expansion needs sampled states.');
end
end


