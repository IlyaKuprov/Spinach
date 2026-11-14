% Evaluate CCE state/cluster tasks in deterministic indexed batches. Syntax:
%
%     raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)
%
% raw is a cell per input context, containing the empty-cluster record
% followed by cluster records in input order. Tasks execute in ascending
% configuration/state/cluster ID order; no reduction occurs on workers.
% parallel.mode is serial (default), clusters or states. Parallel modes
% require a pool already started by create using explicit sys.parallel.
% parallel.batch_size and max_tasks are optional positive integer caps;
% batch_size defaults to the worker count. numerics.max_memory_bytes
% bounds resident raw data plus simultaneous per-worker dense storage.
% checkpoint.file/resume enable per-batch atomic restart. No GPU path.
%
% Inputs:
%   model, clusters, ensemble, schedule, parameters, evaluator - prepared
%   objects and evaluator function handle. Optional parallel.mode,
%   parallel.batch_size, parallel.max_tasks, and checkpoint.file/resume;
%   serial batch_size defaults to one and parallel batch_size to the pool
%   size.
%
% Outputs:
%   raw - cell per input context, each containing the baseline record then
%   cluster records in input order.
%
% See also cce_checkpoint, cce_cluster.
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)

% Validate immutable inputs and choose one outer parallel layer
grumble(model,clusters,ensemble,schedule,parameters);

% Default to the closed-system cluster evaluator
if nargin<6, evaluator=@cce_cluster; end
mode='serial'; batch_size=1; max_tasks=Inf; workers=1;
if isfield(parameters,'parallel')
    if isfield(parameters.parallel,'mode'), mode=parameters.parallel.mode; end
    if isfield(parameters.parallel,'max_tasks'), max_tasks=parameters.parallel.max_tasks; end
end
if ~strcmp(mode,'serial')
    if ~license('test','Distrib_Computing_Toolbox')
        error('cce_schedule:toolbox','Parallel Computing Toolbox is unavailable.');
    end
    pool=gcp('nocreate');
    if isempty(pool)||~isfield(model.spin_system.sys,'parallel')||...
       pool.NumWorkers~=model.spin_system.sys.parallel{2}
        error('cce_schedule:pool','Use create with explicit sys.parallel before scheduling.');
    end
    workers=pool.NumWorkers; batch_size=workers;
end
if isfield(parameters,'parallel')&&isfield(parameters.parallel,'batch_size')
    batch_size=parameters.parallel.batch_size;
end

% Bound task count before allocating the ledger or any cluster operator
ns=numel(ensemble.contexts); nc=numel(clusters)+1; count=ns*nc;
if count>max_tasks, error('cce_schedule:tasks','Task count exceeds parallel.max_tasks.'); end
nt=numel(schedule.times); dc=prod(model.multiplicities(model.central_ids));
components=1;
if strcmp(parameters.quantity,'density'), components=dc^2; end
metadata=whos('model','clusters','ensemble','schedule','parameters');
resident=3*count*(nt*(25*components+16*dc^2)+1024)+sum([metadata.bytes]);

% Include the retained first and second moments of sampled noise records
if strcmp(parameters.method,'noise'), resident=resident+3*count*(32*nt+16); end
peak=0;
for k=0:numel(clusters)
    sites=[];
    if k>0, sites=clusters(k).sites; end
    db=prod(model.multiplicities(sites)); dimension=dc*db;
    if ismember(parameters.method,{'cce','mecce','noise'}), dimension=db; end
    if ismember(parameters.method,{'mecce','megcce'}), dimension=dimension^2; end
    workspace=16*(24+numel(schedule.events))*dimension^2;
    if isfield(schedule,'controls')&&schedule.controls

        % Reserve every retained segment on each simultaneously active worker
        segments=nt*(2+2*numel(schedule.events));
        workspace=workspace+segments*(48*dimension^2+24*numel(model.bath_ids)+8);
    end
    peak=max(peak,workspace);
end
if resident+min([workers batch_size count])*peak>parameters.numerics.max_memory_bytes
    error('cce_schedule:memory','Raw ledger and simultaneous dense storage exceed memory budget.');
end

% Build stable task keys independently of input order and worker rank
ids=zeros(count,3); input_index=zeros(count,2); position=0;
for s=1:ns
    for k=0:numel(clusters)
        position=position+1; id=0;
        if k>0, id=double(clusters(k).id); end
        ids(position,:)=[ensemble.contexts(s).configuration_id ensemble.contexts(s).state_id id];
        input_index(position,:)=[s k];
    end
end
[ids,permutation]=sortrows(ids,[1 2 3]); input_index=input_index(permutation,:);
if size(unique(ids,'rows'),1)~=count
    error('cce_schedule:identity','Configuration, state and cluster IDs must be unique.');
end

% Exclude execution choices from physical restart identity
physical=parameters;
for field={'parallel','checkpoint','retain_raw','convergence'}
    if isfield(physical,field{1}), physical=rmfield(physical,field{1}); end
end

% Hash every executable CCE source so stale code cannot reuse a checkpoint
sources=dir(fullfile(fileparts(mfilename('fullpath')),'*.m'));
code=cell(1,numel(sources));
for k=1:numel(sources), code{k}=fileread(fullfile(sources(k).folder,sources(k).name)); end
identity=struct('model',md5_hash({model.spin_system.comp,model.spin_system.inter}),...
                'parameters',md5_hash(physical),'clusters',md5_hash(clusters),...
                'schedule',md5_hash(schedule),'samples',md5_hash(ensemble),...
                'code',md5_hash(code));
payload=struct('parameters',physical,'ensemble',ensemble,'clusters',clusters,'schedule',schedule);
state=struct('version',1,'identity',identity,'payload',payload,'task_ids',ids,...
             'completed',false(count,1),'records',{cell(count,1)});
filename=''; resume=false;
if isfield(parameters,'checkpoint')
    if isfield(parameters.checkpoint,'file'), filename=parameters.checkpoint.file; end
    if isfield(parameters.checkpoint,'resume'), resume=parameters.checkpoint.resume; end
end
if resume
    state=cce_checkpoint(filename,'read',state);
elseif ~isempty(filename)&&isfile(filename)
    error('cce_schedule:exists','Checkpoint exists; explicitly request resume or use another filename.');
end

% Evaluate only missing tasks and persist each bounded completed batch
missing=find(~state.completed);
for first=1:batch_size:numel(missing)
    indices=missing(first:min(first+batch_size-1,numel(missing)));
    records=cell(numel(indices),1); selected=input_index(indices,:);
    task_clusters=cell(numel(indices),1); task_contexts=cell(numel(indices),1);
    for j=1:numel(indices)
        task_clusters{j}=struct('id',uint64(0),'sites',[]);
        if selected(j,2)>0, task_clusters{j}=clusters(selected(j,2)); end
        task_contexts{j}=ensemble.contexts(selected(j,1));
    end
    if strcmp(mode,'serial')
        for j=1:numel(indices)
            records{j}=evaluator(model,task_clusters{j},task_contexts{j},schedule,parameters);
        end
    else
        parfor (j=1:numel(indices),workers)
            records{j}=evaluator(model,task_clusters{j},task_contexts{j},schedule,parameters); %#ok<PFBNS>
        end
    end
    state.records(indices)=records; state.completed(indices)=true;
    if ~isempty(filename), state=cce_checkpoint(filename,'write',state); end
end

% Restore context/cluster input order for the client's fixed-order reduction
raw=cell(1,ns);
for s=1:ns
    indices=find(input_index(:,1)==s);
    [~,order]=sort(input_index(indices,2));
    raw{s}=[state.records{indices(order)}];
end

end

% Reject unsupported execution fields and malformed scheduling inputs
function grumble(model,clusters,ensemble,schedule,parameters)
if ~isstruct(model)||~all(isfield(model,{'spin_system','multiplicities','central_ids'}))||...
   ~isstruct(clusters)||(~isempty(clusters)&&~all(isfield(clusters,{'id','sites'})))||...
   ~isstruct(ensemble)||~isfield(ensemble,'contexts')||isempty(ensemble.contexts)||...
   ~isstruct(schedule)||~all(isfield(schedule,{'times','events'}))||...
   ~isstruct(parameters)||~all(isfield(parameters,{'numerics','quantity','method'}))||...
   ~isfield(parameters.numerics,'max_memory_bytes')||...
   ~isscalar(parameters.numerics.max_memory_bytes)||...
   ~isfinite(parameters.numerics.max_memory_bytes)||parameters.numerics.max_memory_bytes<=0
    error('cce_schedule:input','Prepared model, clusters, contexts, schedule and memory budget required.');
end
if isfield(parameters,'parallel')
    parallel=parameters.parallel;
    if ~isstruct(parallel)||~isscalar(parallel)||...
       ~all(ismember(fieldnames(parallel),{'mode','batch_size','max_tasks'}))
        error('cce_schedule:input','Unknown parallel field.');
    end
    if isfield(parallel,'mode')&&~ismember(parallel.mode,{'serial','clusters','states'})
        error('cce_schedule:input','parallel.mode must be serial, clusters or states.');
    end
    for field={'batch_size','max_tasks'}
        if isfield(parallel,field{1})
            value=parallel.(field{1});
            if ~isnumeric(value)||~isscalar(value)||~isfinite(value)||value<1||value~=fix(value)
                error('cce_schedule:input','Parallel caps must be finite positive integers.');
            end
        end
    end
end
if isfield(parameters,'checkpoint')
    checkpoint=parameters.checkpoint;
    if ~isstruct(checkpoint)||~isscalar(checkpoint)||...
       ~all(ismember(fieldnames(checkpoint),{'file','resume'}))||...
       (isfield(checkpoint,'file')&&(~ischar(checkpoint.file)||...
        (~isempty(checkpoint.file)&&~isrow(checkpoint.file))))||...
       (isfield(checkpoint,'resume')&&(~islogical(checkpoint.resume)||~isscalar(checkpoint.resume)))
        error('cce_schedule:input','Invalid checkpoint.file/resume.');
    end
    if isfield(checkpoint,'resume')&&checkpoint.resume&&...
       (~isfield(checkpoint,'file')||isempty(checkpoint.file))
        error('cce_schedule:input','Resume requires a checkpoint filename.');
    end
end
end


