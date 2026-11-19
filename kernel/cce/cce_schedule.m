% Evaluate CCE state/cluster tasks in deterministic indexed batches. Syntax:
%
%     raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)
%
% raw is a cell per input spin_config, containing the empty-cluster record
% followed by cluster records in input order. Tasks execute in ascending
% configuration/state/cluster ID order; no reduction occurs on workers.
% parallel.mode is serial (default), clusters, or states. Parallel modes
% require a pool already started by create using explicit sys.parallel.
% parallel.batch_size and max_tasks are optional positive integer caps;
% batch_size defaults to the worker count. numerics.max_memory_bytes
% bounds resident raw data plus simultaneous per-worker dense storage.
% checkpoint.file/resume enable per-batch atomic restart. No GPU path.
% All requests place these controls under parameters.cce; max_tasks caps
% native batch concurrency. Their memory preflight includes retained master
% generators and full physical construction before projection.
%
% Parameters:
%
%   model, clusters, ensemble, schedule, parameters, evaluator - prepared
%   objects and evaluator function handle. Optional parallel.mode,
%   parallel.batch_size, parallel.max_tasks, and checkpoint.file/resume;
%   serial batch_size defaults to one and parallel batch_size to the pool
%   size.
%
% Outputs:
%
%   raw - cell per input spin_config, each containing the baseline record then
%   cluster records in input order.
%
% See also cce_checkpoint, cce_capabilities, clexpn.
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)

% Validate the registered experiment and resolved native scheduling inputs
grumble(model,clusters,ensemble,schedule,parameters,evaluator);

% Build deterministic tasks from resolved native controls
options=parameters.cce;
ns=numel(ensemble.spin_config);
nc=numel(clusters)+1;
mode='serial';
workers=1;
batch_size=1;
if isfield(options,'parallel')
    mode=options.parallel.mode;
end
if ~strcmp(mode,'serial')
    pool=gcp('nocreate');
    if isempty(pool)||~isfield(model.spin_system.sys,'parallel')||...
       pool.NumWorkers~=model.spin_system.sys.parallel{2}
        error('cce_schedule:pool','Parallel scheduling requires an explicitly configured pool.');
    end
    workers=pool.NumWorkers;
    batch_size=workers;
end
if isfield(options,'parallel')&&isfield(options.parallel,'batch_size')
    batch_size=options.parallel.batch_size;
end
if isfield(options,'parallel')&&isfield(options.parallel,'max_tasks')
    limit=options.parallel.max_tasks;
    batch_size=min(batch_size,limit);
end
components=1;
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    components=prod(model.multiplicities(model.observed_ids))^2;
elseif isfield(parameters,'coil')&&iscell(parameters.coil)
    components=numel(parameters.coil);
end
nt=numel(parameters.times);
dc=prod(model.multiplicities(model.observed_ids));
metadata=whos('model','clusters','ensemble','schedule','parameters');
resident=3*ns*nc*(nt*(25*components+16*dc^2)+1024)+sum([metadata.bytes]);
noise=isequal(schedule.pulse_sequence,@noise_correlation);
master=strcmp(model.formalism,'zeeman-liouv')||...
       (isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps));
if noise
    resident=resident+3*ns*nc*(32*nt+16);
end

% Count retained experiment-owned control intervals before dispatch
segment_count=0;
if isequal(schedule.pulse_sequence,@finite_drive)||isequal(schedule.pulse_sequence,@refocusing_train)
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
    ids=[model.observed_ids sites];
    db=prod(model.multiplicities(sites));
    dimension=dc*db;
    copies=12+3*numel(ids);
    if ~model.direct_supported||strcmp(model.backend,'stock')
        ranks=max([2 cellfun(@numel,model.spin_system.inter.giant.coeff).']);
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
    workspace=workspace+16*24*propagated^2+segment_count*(48*dimension^2+24*numel(model.bath_ids)+8);
    peak=max(peak,workspace);
end
if resident+min([workers batch_size ns*nc])*peak>options.numerics.max_memory_bytes
    error('cce_schedule:memory','Raw ledger and simultaneous dense storage exceed memory budget.');
end
model.clusters=clusters;

% Hash physical and experiment inputs, not checkpoint filenames or resume switches
payload=struct('model',model,'clusters',clusters,'ensemble',ensemble,...
               'parameters',parameters,'experiment',func2str(schedule.pulse_sequence),...
               'contract_version',schedule.contract.version);
for field={'checkpoint','parallel','retain_raw'}
    if isfield(payload.parameters.cce,field{1})
        payload.parameters.cce=rmfield(payload.parameters.cce,field{1});
    end
    if isfield(payload.model.parameters.cce,field{1})
        payload.model.parameters.cce=rmfield(payload.model.parameters.cce,field{1});
    end
end

% The model is hashed from its full physical content, not its derived cache key
payload.model=rmfield(payload.model,'source_hash');
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
identity=struct();
if isfield(options,'checkpoint')

    % Available memory is a dispatch limit, not part of the physical calculation
    physical=payload;
    physical.parameters.cce.numerics=rmfield(physical.parameters.cce.numerics,'max_memory_bytes');
    physical.model.parameters.cce.numerics=rmfield(physical.model.parameters.cce.numerics,'max_memory_bytes');
    physical.model.numerics=rmfield(physical.model.numerics,'max_memory_bytes');
    identity=struct('model',md5_hash(physical.model),'parameters',md5_hash(physical.parameters),...
                    'clusters',md5_hash(clusters),'schedule',md5_hash(payload.experiment),...
                    'samples',md5_hash(ensemble),'code',payload.contract_version);

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
for b=1:batch_size:numel(pending)
    selected=pending(b:min(b+batch_size-1,numel(pending)));
    records=cell(numel(selected),1);
    tasks=indices(selected,:);
    if strcmp(mode,'serial')
        for j=1:numel(selected)
            records{j}=native_task(model,clusters,ensemble,schedule,parameters,evaluator,tasks(j,:));
        end
    else
        parfor (j=1:numel(selected),workers)
            records{j}=native_task(model,clusters,ensemble,schedule,parameters,evaluator,tasks(j,:));
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
        position=find(indices(:,1)==n&indices(:,2)==k,1);
        records{k}=state.records{position};
    end
    if strcmp(options.bath_state.mode,'interlaced')
        records=records(~cellfun(@(item)strcmp(item.kind,'not_evaluated'),records));
    end
    raw{n}=[records{:}];
end

end

% Invoke exactly the registered evaluator on ordinary and conditioned paths
function raw=native_task(model,clusters,ensemble,schedule,parameters,evaluator,index)
cluster=struct('id',0,'sites',[]);
if index(2)>1
    cluster=clusters(index(2)-1);
end
if strcmp(parameters.cce.bath_state.mode,'interlaced')
    parameters.cce.interlaced.reference_state=index(1);
    complete=any(arrayfun(@(item)isequal(sort(item.sites),sort(model.bath_ids)),clusters));
    if complete&&~isequal(sort(cluster.sites),sort(model.bath_ids))
        channels=1;
        if isfield(parameters,'coil')&&iscell(parameters.coil)
            channels=numel(parameters.coil);
        end
        raw=struct('cluster_id',cluster.id,'value',NaN(numel(parameters.times),channels),...
                   'abs_error_bound',Inf(numel(parameters.times),channels),...
                   'valid',false(numel(parameters.times),channels),...
                   'weight',ensemble.spin_config(index(1)).weight,...
                   'condition_sites',{{}},'assignment_count',0,'reference',[],...
                   'reference_error',[],'full_bath',false,'kind','not_evaluated');
        if strcmp(model.formalism,'zeeman-liouv')||...
           (isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps))
            raw.generator=[];
        end
    else
        raw=cce_interlace(model,cluster,ensemble,schedule,parameters,evaluator);
    end
else
    raw=evaluator(model,cluster,ensemble.spin_config(index(1)),schedule,parameters);
end
end

% Validate producer-owned records and the required explicit callback
function grumble(model,clusters,ensemble,schedule,parameters,evaluator)
if ~isstruct(model)||~isscalar(model)||~all(isfield(model,{'spin_system','multiplicities','observed_ids'}))||...
   ~isstruct(clusters)||(~isempty(clusters)&&~all(isfield(clusters,{'id','sites'})))||...
   ~isstruct(ensemble)||~isfield(ensemble,'spin_config')||isempty(ensemble.spin_config)||...
   ~isstruct(schedule)||~all(isfield(schedule,{'pulse_sequence','contract'}))||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'cce')||...
   ~isa(evaluator,'function_handle')
    error('cce_schedule:input','Native model, ensemble, schedule, parameters and evaluator are required.');
end
mode='serial';
if isfield(parameters.cce,'parallel')
    mode=parameters.cce.parallel.mode;
end
if ~ismember(mode,{'serial','clusters','states'})
    error('cce_schedule:parallel','Choose serial, clusters or states parallelism.');
end
if isfield(parameters.cce,'parallel')&&isfield(parameters.cce.parallel,'batch_size')
    batch_size=parameters.cce.parallel.batch_size;
    if ~isscalar(batch_size)||~isfinite(batch_size)||batch_size<1||batch_size~=fix(batch_size)
        error('cce_schedule:batch','Batch size must be a finite positive integer.');
    end
end
if isfield(parameters.cce,'parallel')&&isfield(parameters.cce.parallel,'max_tasks')
    limit=parameters.cce.parallel.max_tasks;
    if ~isnumeric(limit)||~isreal(limit)||~isscalar(limit)||...
       ~isfinite(limit)||limit<1||limit~=fix(limit)
        error('parallel.max_tasks must be a positive integer.');
    end
end
if strcmp(parameters.cce.bath_state.mode,'interlaced')&&...
   ~isfield(parameters.cce,'interlaced')
    error('cce_schedule:interlaced','Interlaced expansion needs sampled states and an enumeration budget.');
end
end


