% Cluster correlation expansion for a prepared Spinach spin system. Syntax:
%
%                       result=cce(spin_system,parameters)
%
% Required fields: observed_spins, bath_spins, method (cce/gcce/mecce/megcce/noise), expansion
% (standard), quantity (signal/density/observable/autocorrelation), times (whole-sequence
% seconds), order, observed, assumptions, orientation, interaction_cutoff,
% graph, bath_state, meanfield, sequence, numerics.atol/rtol/max_memory_bytes.
% cce_prepare, cce_graph, cce_states and cce_sequence document their schemas.
% normalized=true divides gcce observables by the supplied initial state;
% projected manifold signal is intrinsically normalized. Density is not.
% reducer is powers (default) or direct. Sampled estimators explicitly set
% averaging='after_expansion' or 'internal_then_expansion'.
%
% Outputs include complex value/signal, per-order estimates, validity,
% reasons, propagated error bounds, raw clusters and exact-full-bath bypass
% metadata. Density outputs additionally report physicality without repair.
% Extensions fail with cce:capability until their implementations are present;
% no optional control, checkpoint or dissipative physics is silently ignored.
% ME jumps use ordered physical support, dimensionless operator and s^-1 rate
% fields (cce_liouvillian). Noise requires explicit ordering, connected flag
% and Hermitian support/operator terms (cce_noise); its reducer is additive.
%
% Inputs:
%   spin_system - canonical structure returned by create.
%   parameters - explicit physical, state, graph, sequence and numerical
%   choices; see the CCE recipe.
%
% Outputs:
%   result - reduction values/masks/error/by-order fields plus signal,
%   times, raw, clusters, graph, schedule, partition, ensemble, model_hash,
%   input_parameters, averaging, method, exact_reference, elapsed_seconds,
%   convergence, budget_satisfied, sampling_standard_error and diagnostics;
%   density is empty except for density requests. value is time-by-observable;
%   density elements use column-major order.
%
% Multiplicative convention: F(empty)=1; F(C)=L(C)/L(empty),
% K(C)=F(C)/product(K(S), nonempty proper S subset C), and
% L_approx=L(empty)*product(K(C), retained nonempty C).
% For noise, K(C)=L(C)-L(empty)-sum(K(S), proper nonempty S),
% and L_approx=L(empty)+sum(K(C)). Masks remain explicit;
% unresolved denominators are never replaced by epsilon.
%
% See also cce_prepare, cce_cluster, cce_reduce, cce_diagnostics.
% Reference: https://doi.org/10.1103/PhysRevB.78.085315
% Reference: https://doi.org/10.1103/PRXQuantum.2.010311
%
% ilya.kuprov@weizmann.ac.il

function result=cce(spin_system,parameters)

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Reject conditioning controls when the driver is not performing interlacing
if isfield(parameters,'interlaced')&&isfield(parameters,'bath_state')&&...
   ~strcmp(parameters.bath_state.mode,'interlaced')
    error('cce:capability','parameters.interlaced requires bath_state.mode=interlaced.');
end

% Convergence scans are explicit helper operations, not driver flags
if isfield(parameters,'convergence')
    error('cce:capability','parameters.convergence requires explicit cce_converge execution.');
end

% Validate the core driver contract before allocating retained results
grumble(spin_system,parameters);
started=tic; model=cce_prepare(spin_system,parameters);

% Dispatch open-system and additive-noise cluster evaluators
evaluator=@cce_cluster;
if ismember(parameters.method,{'mecce','megcce'}), evaluator=@cce_master; end
if strcmp(parameters.method,'noise'), evaluator=@cce_noise; end
state_parameters=parameters;
if ismember(parameters.method,{'mecce','noise'})
    state_parameters.method='cce';
    if isfield(state_parameters,'jumps'), state_parameters=rmfield(state_parameters,'jumps'); end
    if isfield(state_parameters.numerics,'propagator')
        state_parameters.numerics=rmfield(state_parameters.numerics,'propagator');
    end
end
if strcmp(parameters.method,'noise')
    state_parameters=rmfield(state_parameters,'noise');
    if isfield(state_parameters,'spectrum'), state_parameters=rmfield(state_parameters,'spectrum'); end
end
ensemble=cce_states(model,state_parameters); schedule=cce_sequence(model,parameters);
graph=cce_graph(model,parameters);

% Freeze physical group partitions before evaluating any production sample
partition=[];
if ~strcmp(parameters.expansion,'standard')
    partition=cce_partition(model,graph,parameters);
    if ~partition.complete, error('cce:partition_incomplete','Partition validation did not converge.'); end
    parameters.groups=partition.groups;
    if strcmp(parameters.expansion,'disjoint'), parameters.order=1; end
    if strcmp(parameters.expansion,'disjoint_inter'), parameters.order=2; end
end
clusters=cce_clusters(model,graph,parameters);
full_index=find(arrayfun(@(c) isequal(sort(c.sites),sort(model.bath_ids)),clusters),1);
if isempty(model.bath_ids), full_index=0; end

% Apply W12 conditional averages before irreducible subset reduction
if strcmp(parameters.bath_state.mode,'interlaced')
    if ~strcmp(parameters.quantity,'signal')
        error('cce:interlaced','Interlaced multiplication requires scalar signal.');
    end
    empty=struct('id',uint64(0),'sites',[]); ns=numel(ensemble.spin_configs);
    nt=numel(parameters.times); no=parameters.order+1;

    % Bound retained sample/order values, uncertainties and validity masks
    dc=prod(model.multiplicities(model.observed_ids));
    bytes=3*nt*ns*(25+25*no)+...
          3*ns*(numel(clusters)+1)*(nt*(25+16*dc^2)+1024);
    if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
        error('cce:memory',...
              'Interlaced sample and order summaries exceed numerics.max_memory_bytes.');
    end
    estimates=zeros(nt,ns); by_order=zeros(nt,1,no,ns);
    order_error=zeros(size(by_order)); order_valid=true(size(by_order));
    valid=true(size(estimates)); bounds=zeros(size(estimates));

    % Retained conditioning records scale with both samples and cluster count
    cylinders={clusters.sites};
    if isfield(parameters.interlaced,'condition_sites')
        cylinders=parameters.interlaced.condition_sites;
    end
    if isempty(cylinders), cylinders={[]}; end
    conditioning=0;
    for sites=cylinders
        conditioning=conditioning+prod(model.multiplicities(sites{1}));
    end
    conditioning=min(conditioning,parameters.interlaced.max_states);
    count=ns*(numel(clusters)+1);
    if ~isempty(full_index), count=ns; end
    retained=16*conditioning*(numel(model.bath_ids)+1)*count;
    if ~isfinite(retained)||bytes+retained>parameters.numerics.max_memory_bytes
        error('cce:memory',...
              'Retained conditioning records exceed numerics.max_memory_bytes.');
    end
    records=cell(ns,numel(clusters)+1);
    for s=1:ns
        conditioned=parameters; conditioned.interlaced.reference_state=s;
        if ~isempty(full_index)

            % Full-bath conditioning is an exact raw average without denominators
            conditioned.interlaced.raw_only=true;
            if full_index==0, cluster=empty; else, cluster=clusters(full_index); end
            item=cce_interlace(model,cluster,ensemble,schedule,conditioned);
            estimates(:,s)=item.value; bounds(:,s)=item.abs_error_bound;
            valid(:,s)=item.valid; records{s,full_index+1}=item;
            by_order(:,:,1:end-1,s)=NaN; order_error(:,:,1:end-1,s)=Inf;
            order_valid(:,:,1:end-1,s)=false;
            by_order(:,:,end,s)=item.value; order_error(:,:,end,s)=item.abs_error_bound;
            order_valid(:,:,end,s)=item.valid; continue;
        end
        base=cce_interlace(model,empty,ensemble,schedule,conditioned);
        product=base.value; uncertainty=base.abs_error_bound; valid(:,s)=base.valid;
        records{s,1}=base; by_order(:,:,1,s)=product;
        order_error(:,:,1,s)=uncertainty; order_valid(:,:,1,s)=valid(:,s);
        for order=1:parameters.order
            for j=find([clusters.order]==order)
                item=cce_interlace(model,clusters(j),ensemble,schedule,conditioned);
                uncertainty=abs(product).*item.factor_abs_error_bound+...
                    abs(item.factor).*uncertainty+uncertainty.*item.factor_abs_error_bound;
                product=product.*item.factor; valid(:,s)=valid(:,s)&item.factor_valid;
                records{s,j+1}=item;
            end
            by_order(:,:,order+1,s)=product; order_error(:,:,order+1,s)=uncertainty;
            order_valid(:,:,order+1,s)=valid(:,s);
        end
        estimates(:,s)=product; bounds(:,s)=uncertainty;
    end

    % Aggregate every computed prefix and preserve unavailable lower exact orders
    weights=reshape(ensemble.weights,1,1,1,[]);
    value=estimates*ensemble.weights(:); uncertainty=bounds*ensemble.weights(:);
    prefixes=sum(by_order.*weights,4); prefix_error=sum(order_error.*weights,4);
    prefix_valid=all(order_valid,4); reasons=repmat({'ok'},size(prefixes));
    reasons(~prefix_valid)={'unresolved_interlaced_factor'};
    if ~isempty(full_index), reasons(:,:,1:end-1)={'not_evaluated'}; end
    mask=all(valid,2)&isfinite(value)&isfinite(uncertainty);
    result=struct('value',value,'valid',mask,'abs_error_bound',uncertainty,...
                  'reason',{reasons(:,:,end)},'by_order',prefixes,...
                  'by_order_error',prefix_error,'by_order_valid',prefix_valid,...
                  'by_order_reason',{reasons},'orders',0:parameters.order,'reducer','interlaced');
    result.raw=records(~cellfun(@isempty,records));
    result.sampling_standard_error=NaN(size(value));
    if ns>1, result.sampling_standard_error=std(estimates,0,2)/sqrt(ns); end
    if ~isempty(full_index), result.sampling_standard_error(:)=0; end
    averaging='conditional_interlaced';
else

% Propagate each globally assigned bath state before ensemble reduction
ns=numel(ensemble.spin_configs); reductions=cell(1,ns);
all_raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator);

% Preserve the covariance of the sampled mixture rather than conditional covariances
if strcmp(parameters.method,'noise')&&parameters.noise.connected&&ns>1

    % Average cluster moments before taking the connected covariance
    for k=1:numel(all_raw{1})
        mean_time=zeros(size(all_raw{1}(k).value)); mean_zero=0;
        for s=1:ns
            item=all_raw{s}(k); weight=ensemble.weights(s);
            mean_time=mean_time+weight*item.noise_mean;
            mean_zero=mean_zero+weight*item.noise_mean_zero;
        end

        % Retain linear sample contributions for either additive estimator
        for s=1:ns
            all_raw{s}(k).value=all_raw{s}(k).noise_moment-...
                                all_raw{s}(k).noise_mean*mean_zero-...
                                mean_time*(all_raw{s}(k).noise_mean_zero-mean_zero);
        end
    end
end

% Distinguish average-of-expansions from an explicitly requested internal average
averaging='after_expansion';
if isfield(parameters,'averaging'), averaging=parameters.averaging; end
if strcmp(averaging,'internal_then_expansion')&&ns>1
    records=all_raw{1};
    for k=1:numel(records)
        records(k).value(:)=0; records(k).abs_error_bound(:)=0;
        for s=1:ns
            records(k).value=records(k).value+ensemble.weights(s)*all_raw{s}(k).value;
            records(k).abs_error_bound=records(k).abs_error_bound+...
                                      ensemble.weights(s)*all_raw{s}(k).abs_error_bound;
            records(k).valid=records(k).valid&all_raw{s}(k).valid;
        end
    end
    reduction_inputs={records}; weights=1;
else
    reduction_inputs=all_raw; weights=ensemble.weights;
end

% Use exact full-bath raw values even when intermediate baselines vanish
for s=1:numel(reduction_inputs)
    records=reduction_inputs{s};
    reduction_parameters=parameters;
    if ~isempty(full_index), reduction_parameters.numerics.zero_policy='mask'; end
    reduced=cce_reduce(records,clusters,reduction_parameters);
    if ~isempty(full_index)
        exact=records(full_index+1);
        reduced.value=exact.value; reduced.valid=exact.valid;
        reduced.abs_error_bound=exact.abs_error_bound;
        reduced.reason=repmat({'ok'},size(exact.value));
        reduced.by_order(:,:,end)=exact.value;
        reduced.by_order_valid(:,:,end)=exact.valid;
        reduced.by_order_error(:,:,end)=exact.abs_error_bound;
        reduced.by_order_reason(:,:,end)=repmat({'ok'},size(exact.value));
    end
    reductions{s}=reduced;
end

% Combine completed samples deterministically without masking failed samples
result=reductions{1}; result.value(:)=0; result.abs_error_bound(:)=0;
result.by_order(:)=0; result.by_order_error(:)=0;
for s=1:numel(reduction_inputs)
    result.value=result.value+weights(s)*reductions{s}.value;
    result.abs_error_bound=result.abs_error_bound+weights(s)*reductions{s}.abs_error_bound;
    result.by_order=result.by_order+weights(s)*reductions{s}.by_order;
    result.by_order_error=result.by_order_error+weights(s)*reductions{s}.by_order_error;
    result.by_order_valid=result.by_order_valid&reductions{s}.by_order_valid;
    failed_order=~reductions{s}.by_order_valid;
    result.by_order_reason(failed_order)=reductions{s}.by_order_reason(failed_order);
    result.valid=result.valid&reductions{s}.valid;
    failed=~reductions{s}.valid; result.reason(failed)=reductions{s}.reason(failed);
end

% Account for deterministic weighted summation without hiding sample uncertainty
if numel(reduction_inputs)>1
    operations=4*numel(reduction_inputs);
    rounding=operations*eps/(1-operations*eps); magnitude=zeros(size(result.value));
    for s=1:numel(reduction_inputs)
        magnitude=magnitude+weights(s)*abs(reductions{s}.value);
    end
    result.abs_error_bound=result.abs_error_bound+rounding*magnitude;
end

% Expose estimator uncertainty separately from propagation error bounds
result.sampling_standard_error=zeros(size(result.value));
if strcmp(ensemble.mode,'sampled')
    result.sampling_standard_error(:)=NaN; values=[];
    if ns>1&&strcmp(averaging,'after_expansion')
        values=cellfun(@(item)item.value,reductions,'UniformOutput',false);
        values=cat(3,values{:});
    elseif ns>1&&~isempty(full_index)&&...
           ~(strcmp(parameters.method,'noise')&&parameters.noise.connected)
        values=cellfun(@(records)records(full_index+1).value,all_raw,'UniformOutput',false);
        values=cat(3,values{:});
    end
    if ~isempty(values)
        result.sampling_standard_error=sqrt(sum(abs(values-mean(values,3)).^2,3)/(ns*(ns-1)));
    end
end
result.raw=all_raw;
end

% Finalize every estimator through the same public result and diagnostic contract
result.response=result.value; result.times=parameters.times;
result.clusters=clusters; result.graph=graph; result.schedule=schedule;
result.partition=partition;
result.ensemble=ensemble; result.model_hash=model.source_hash;
result.input_parameters=parameters; result.averaging=averaging;
result.method=parameters.method; result.exact_reference=~isempty(full_index);
if result.exact_reference, result.method='exact-reference'; end
result.elapsed_seconds=toc(started); result.convergence='not_requested';
result.budget_satisfied=result.valid&...
    result.abs_error_bound<=parameters.numerics.atol+parameters.numerics.rtol*abs(result.value);

% Diagnose a reconstructed observed density without symmetrization or clipping
result.density=[];
if strcmp(parameters.quantity,'density')
    d=ensemble.observed_dimension; nt=numel(parameters.times);
    result.density=reshape(result.value.',d,d,nt);
end
result.diagnostics=cce_diagnostics(model,result.raw,result,parameters);

% Apply explicitly requested postprocessing to the computed result
if isfield(parameters,'spectrum')
    result.spectrum=cce_spectrum(result.times,result.value,parameters);
end
if isfield(parameters,'analysis'), result.analysis=cce_analyse(result,parameters); end

% Honor explicit raw suppression after diagnostics and analysis consume it
if isfield(parameters,'retain_raw')&&~parameters.retain_raw, result.raw={}; end

end

% Retained driver-local scheduling pending removal of this entry point
function raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate immutable inputs and choose one outer parallel layer
schedule_grumble(model,clusters,ensemble,schedule,parameters);

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
ns=numel(ensemble.spin_configs); nc=numel(clusters)+1; count=ns*nc;
if count>max_tasks, error('cce_schedule:tasks','Task count exceeds parallel.max_tasks.'); end
nt=numel(schedule.times); dc=prod(model.multiplicities(model.observed_ids));
components=1;
if strcmp(parameters.quantity,'density'), components=dc^2; end
if strcmp(parameters.quantity,'observable')&&iscell(parameters.coil)
    components=numel(parameters.coil);
end
metadata=whos('model','clusters','ensemble','schedule','parameters');
resident=3*count*(nt*(25*components+16*dc^2)+1024)+sum([metadata.bytes]);

% Include the retained first and second moments of sampled noise records
if strcmp(parameters.method,'noise'), resident=resident+3*count*(32*nt+16); end
peak=0;
for k=0:numel(clusters)
    sites=[];
    if k>0, sites=clusters(k).sites; end

    % Reserve unprojected Hamiltonians, spin operators and physical pulse matrices
    db=prod(model.multiplicities(sites)); dimension=dc*db;
    copies=12+3*(numel(model.observed_ids)+numel(sites))+2*numel(schedule.events);
    if ~model.direct_supported||strcmp(model.backend,'stock')

        % Include the spherical components retained by the stock construction route
        ranks=max([2 cellfun(@numel,model.spin_system.inter.giant.coeff).']);
        copies=copies+sum((2*(1:ranks)+1).^2);
    end
    workspace=16*copies*dimension^2;

    % Add the projected Hilbert or Liouville propagation workspace separately
    propagated=dimension;
    if ismember(parameters.method,{'cce','mecce','noise'}), propagated=db; end
    if ismember(parameters.method,{'mecce','megcce'})

        % Retain manifold, spin and included jump operators in every completed record
        ids=[model.observed_ids sites]; jumps=0; labels=0;
        for j=1:numel(parameters.jumps)
            if all(ismember(parameters.jumps(j).support,ids))
                jumps=jumps+1; labels=labels+numel(parameters.jumps(j).support)+1;
            end
        end
        retained=16*(3*numel(ids)*dimension^2+2*(1+jumps)*propagated^2)+...
                 8*(2*numel(ids)+labels);
        resident=resident+3*ns*retained; workspace=workspace+retained;
        propagated=propagated^2;
    end
    workspace=workspace+16*(24+numel(schedule.events))*propagated^2;
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
        ids(position,:)=[ensemble.spin_configs(s).configuration_id ensemble.spin_configs(s).state_id id];
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
identity=struct('model',md5_hash({model.spin_system.comp,model.spin_system.inter,...
                                model.spin_system.tols}),...
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

% Keep scheduler-owned checkpoint options out of direct evaluator calls
if isfield(parameters,'checkpoint'), parameters=rmfield(parameters,'checkpoint'); end

% Evaluate only missing tasks and persist each bounded completed batch
missing=find(~state.completed);
for first=1:batch_size:numel(missing)
    indices=missing(first:min(first+batch_size-1,numel(missing)));
    records=cell(numel(indices),1); selected=input_index(indices,:);
    task_clusters=cell(numel(indices),1); task_configs=cell(numel(indices),1);
    for j=1:numel(indices)
        task_clusters{j}=struct('id',uint64(0),'sites',[]);
        if selected(j,2)>0, task_clusters{j}=clusters(selected(j,2)); end
        task_configs{j}=ensemble.spin_configs(selected(j,1));
    end
    if strcmp(mode,'serial')
        for j=1:numel(indices)
            records{j}=evaluator(model,task_clusters{j},task_configs{j},schedule,parameters);
        end
    else
        parfor (j=1:numel(indices),workers)
            records{j}=evaluator(model,task_clusters{j},task_configs{j},schedule,parameters); %#ok<PFBNS>
        end
    end
    state.records(indices)=records; state.completed(indices)=true;
    if ~isempty(filename), state=cce_checkpoint(filename,'write',state); end
end

% Restore spin_config/cluster input order for the client's fixed-order reduction
raw=cell(1,ns);
for s=1:ns
    indices=find(input_index(:,1)==s);
    [~,order]=sort(input_index(indices,2));
    raw{s}=[state.records{indices(order)}];
end

end

% Reject unsupported execution fields and malformed scheduling inputs
function schedule_grumble(model,clusters,ensemble,schedule,parameters)
if ~isstruct(model)||~all(isfield(model,{'spin_system','multiplicities','observed_ids'}))||...
   ~isstruct(clusters)||(~isempty(clusters)&&~all(isfield(clusters,{'id','sites'})))||...
   ~isstruct(ensemble)||~isfield(ensemble,'spin_configs')||isempty(ensemble.spin_configs)||...
   ~isstruct(schedule)||~all(isfield(schedule,{'times','events'}))||...
   ~isstruct(parameters)||~all(isfield(parameters,{'numerics','quantity','method'}))||...
   ~isfield(parameters.numerics,'max_memory_bytes')||...
   ~isscalar(parameters.numerics.max_memory_bytes)||...
   ~isfinite(parameters.numerics.max_memory_bytes)||parameters.numerics.max_memory_bytes<=0
    error('cce_schedule:input','Prepared model, clusters, spin_configs, schedule and memory budget required.');
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

% Retained driver-local scheduling pending removal of this entry point
function raw=cce_interlace(model,cluster,ensemble,schedule,parameters,evaluator)

% Preserve an explicitly supplied experiment evaluator through every recursion
if nargin<6
    evaluator=@cce_cluster;
    if ismember(parameters.method,{'mecce','megcce'}), evaluator=@cce_master; end
end
if ~isa(evaluator,'function_handle')
    error('cce_interlace:evaluator','evaluator must be a function handle.');
end

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate the explicit conditioning and local product distribution
interlace_grumble(model,cluster,ensemble,parameters);

% Locate included superclusters unless a shared denominator K was supplied
spec=parameters.interlaced; reference=ensemble.spin_configs(spec.reference_state);
if isfield(reference,'assignment')&&...
   ~isequal(reference.assignment,ensemble.assignments(spec.reference_state,:))
    error('cce:spin_config','Interlaced spin_config and ensemble reference assignments must agree.');
end
reference.assignment=ensemble.assignments(spec.reference_state,:);
bases=cell(size(model.bath_ids)); probabilities=bases;
for j=1:numel(model.bath_ids)
    d=model.multiplicities(model.bath_ids(j));
    if ischar(parameters.bath_state.basis), bases{j}=eye(d); else, bases{j}=parameters.bath_state.basis{j}; end
    if ischar(parameters.bath_state.probabilities)
        probabilities{j}=ones(d,1)/d;
    else
        probabilities{j}=parameters.bath_state.probabilities{j}(:);
    end
end

% Validate the supplied reference before replacing any enumerated bath state
cce_validate(parameters,model,cluster.sites,reference);
if ~isequal(reference.bath_ids,model.bath_ids)||...
   numel(reference.assignment)~=numel(model.bath_ids)
    error('cce:spin_config','Interlaced reference bath IDs and assignments must match the model.');
end

% Require each reference density to represent its declared local assignment
for j=1:numel(model.bath_ids)
    dimension=model.multiplicities(model.bath_ids(j)); basis=bases{j};
    index=reference.assignment(j);
    if ~isfinite(index)||index~=fix(index)||index<1||index>dimension
        error('cce:spin_config','Interlaced reference assignments must select local basis states.');
    end
    vector=basis(:,index);
    if ~isequal(size(reference.bath_rho{j}),[dimension dimension])||...
       norm(reference.bath_rho{j}-vector*vector','fro')>dimension^2*eps
        error('cce:spin_config',...
              'Interlaced reference bath_rho must agree with its basis assignment; use explicit_product for mixed states.');
    end
end

% Resolve the cylinder family only after the reference has been validated
if isfield(spec,'condition_sites')
    cylinders=spec.condition_sites;
else
    family=cce_clusters(model,cce_graph(model,parameters),parameters);

    % Retain the reference-assignment baseline when no clusters were kept
    if isempty(family)&&isempty(cluster.sites)
        cylinders={[]};
    else
        cylinders={family(arrayfun(@(c) all(ismember(cluster.sites,c.sites)),family)).sites};
    end
end
if isempty(cylinders), error('cce_interlace:empty','No conditioning supercluster.'); end

% Bound cylinder masks, enumeration counts and index vectors before allocation
bytes=3*8*numel(cylinders);
if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
    error('cce_interlace:memory',...
          'Conditioning cylinder metadata exceed numerics.max_memory_bytes.');
end

% Remove contained cylinders without changing the conditioning union
keep=true(size(cylinders));
for j=1:numel(cylinders)
    for k=1:numel(cylinders)
        if j~=k&&all(ismember(cylinders{j},cylinders{k}))&&...
           (numel(cylinders{j})<numel(cylinders{k})||j>k)
            keep(j)=false; break;
        end
    end
end
cylinders=cylinders(keep);

% Reject invalid cylinder IDs and bound enumeration before allocating it
counts=zeros(size(cylinders)); dimensions=model.multiplicities(model.bath_ids);
for k=1:numel(cylinders)
    if ~isnumeric(cylinders{k})||~isreal(cylinders{k})||numel(unique(cylinders{k}))~=numel(cylinders{k})||...
       any(~ismember(cylinders{k},model.bath_ids))||~all(ismember(cluster.sites,cylinders{k}))
        error('cce_interlace:cylinder','Invalid conditioning supercluster.');
    end
    counts(k)=prod(dimensions(ismember(model.bath_ids,cylinders{k})));
end
if sum(counts)>spec.max_states
    error('cce_interlace:budget','Cylinder enumeration exceeds interlaced.max_states.');
end

% Bound assignments, a cylinder block, sorting copies and index/weight vectors
rows=sum(counts); columns=numel(model.bath_ids);
bytes=8*rows*(4*columns+5);
if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
    error('cce_interlace:memory',...
          'Conditioning assignments and sorting workspace exceed numerics.max_memory_bytes.');
end

% Enumerate each cylinder around the same global assignment
assignments=zeros(rows,columns); cursor=0;
for k=1:numel(cylinders)
    positions=find(ismember(model.bath_ids,cylinders{k}));
    block=repmat(reference.assignment,counts(k),1); radix=1;
    for j=positions
        block(:,j)=mod(floor((0:counts(k)-1)'/radix),dimensions(j))+1;
        radix=radix*dimensions(j);
    end
    assignments(cursor+(1:counts(k)),:)=block; cursor=cursor+counts(k);
end

% Release the last cylinder before sorting the conditioning union
clear block

% Deduplicate overlapping cylinders before computing physical probabilities
assignments=unique(assignments,'rows','sorted');
logweights=zeros(size(assignments,1),1);
for j=1:numel(model.bath_ids)
    logweights=logweights+log(probabilities{j}(assignments(:,j)));
end
shift=max(logweights);
if ~isfinite(shift), error('cce_interlace:zero_mass','Conditioning event has zero probability.'); end
weights=exp(logweights-shift); logmass=shift+log(sum(weights));
weights=weights/sum(weights); mass=exp(logmass);

% Evaluate each distinct configuration once using the unchanged cluster engine
raw=[];
for k=1:numel(weights)
    if weights(k)==0, continue; end
    spin_config=reference; spin_config.assignment=assignments(k,:);
    for j=1:numel(model.bath_ids)
        vector=bases{j}(:,assignments(k,j));
        spin_config.bath_rho{j}=vector*vector';
        spin=pauli(dimensions(j)); x=spin.x; y=spin.y; z=spin.z; operators={x,y,z};
        for a=1:3
            spin_config.expectations(j,a)=real(trace(spin_config.bath_rho{j}*operators{a}));
        end
    end

    % Preserve the physical evaluator in every recursive conditional average
    item=evaluator(model,cluster,spin_config,schedule,parameters);
    if isempty(raw)
        raw=item; raw.value=zeros(size(item.value));
        raw.abs_error_bound=zeros(size(item.abs_error_bound));
        raw.density=zeros(size(item.density));
    end
    raw.value=raw.value+weights(k)*item.value;
    raw.abs_error_bound=raw.abs_error_bound+weights(k)*item.abs_error_bound;
    raw.density=raw.density+weights(k)*item.density;
    raw.valid=raw.valid&item.valid;

    % Preserve worst constituent residuals rather than only the first assignment
    if isfield(item,'eigen_residual')
        raw.eigen_residual=max(raw.eigen_residual,item.eigen_residual);
    end
    if isfield(item,'orthogonality_residual')
        raw.orthogonality_residual=max(raw.orthogonality_residual,item.orthogonality_residual);
    end
    if isfield(item,'unitarity_residual')
        raw.unitarity_residual=max(raw.unitarity_residual,item.unitarity_residual);
    end
end

% Form the irreducible W12 factor using precisely the same conditioning event
if (~isfield(spec,'raw_only')||~spec.raw_only)&&~isempty(cluster.sites)
    family=cce_clusters(model,cce_graph(model,parameters),parameters);
    subsets=family(arrayfun(@(entry) all(ismember(entry.sites,cluster.sites)),family));
    bytes=8*rows*(columns+1)+3*25*numel(raw.value)*numel(subsets);
    if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
        error('cce_interlace:memory',...
              'Conditional subset factors exceed numerics.max_memory_bytes.');
    end
    factors=cell(1,numel(subsets)); errors=factors; masks=factors;
    conditioned=parameters; conditioned.interlaced.condition_sites=cylinders;
    conditioned.interlaced.raw_only=true;
    empty=cce_interlace(model,struct('id',0,'sites',[]),ensemble,schedule,conditioned,evaluator);
    for j=1:numel(subsets)
        if isequal(subsets(j).sites,cluster.sites)
            item=raw;
        else
            item=cce_interlace(model,subsets(j),ensemble,schedule,conditioned,evaluator);
        end
        denominator=empty.value; uncertainty=empty.abs_error_bound; valid=empty.valid;
        for k=1:j-1
            if all(ismember(subsets(k).sites,subsets(j).sites))
                uncertainty=abs(denominator).*errors{k}+abs(factors{k}).*uncertainty+uncertainty.*errors{k};
                denominator=denominator.*factors{k}; valid=valid&masks{k};
            end
        end
        resolved=valid&item.valid&abs(denominator)>uncertainty+parameters.numerics.atol+...
                 parameters.numerics.rtol*abs(denominator);
        factors{j}=NaN(size(item.value)); errors{j}=Inf(size(item.value)); masks{j}=resolved;
        factors{j}(resolved)=item.value(resolved)./denominator(resolved);
        errors{j}(resolved)=(item.abs_error_bound(resolved)+...
            abs(factors{j}(resolved)).*uncertainty(resolved))./...
            (abs(denominator(resolved))-uncertainty(resolved));
    end
    target=find(arrayfun(@(entry) isequal(sort(entry.sites),sort(cluster.sites)),subsets));
    if numel(target)~=1, error('cce_interlace:family','Requested cluster absent from included family.'); end
    raw.factor=factors{target}; raw.factor_valid=masks{target}&isfinite(raw.factor);
    raw.factor_abs_error_bound=errors{target};
end

% Return the exact conditioning signature for safe reuse by denominator subsets
raw.condition_sites=cylinders; raw.assignments=assignments;
raw.weights=weights; raw.condition_mass=mass; raw.log_condition_mass=logmass;
raw.distribution_signature=struct('assignments',assignments,'weights',weights);

end

function interlace_grumble(model,cluster,ensemble,parameters)
if ~isstruct(model)||~all(isfield(model,{'bath_ids','multiplicities'}))||...
   ~isstruct(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   ~isstruct(ensemble)||~all(isfield(ensemble,{'spin_configs','assignments'}))||...
   ~isstruct(parameters)||~isfield(parameters,'interlaced')||...
   ~isfield(parameters,'bath_state')||~all(isfield(parameters.bath_state,{'basis','probabilities'}))
    error('cce_interlace:input','Canonical model, cluster, ensemble and interlaced options required.');
end
spec=parameters.interlaced;
if ~isstruct(spec)||~all(isfield(spec,{'reference_state','max_states'}))||...
   ~isscalar(spec.reference_state)||~isreal(spec.reference_state)||...
   ~isfinite(spec.reference_state)||spec.reference_state<1||...
   spec.reference_state~=fix(spec.reference_state)||spec.reference_state>numel(ensemble.spin_configs)||...
   ~isscalar(spec.max_states)||~isreal(spec.max_states)||~isfinite(spec.max_states)||...
   spec.max_states<1||spec.max_states~=fix(spec.max_states)
    error('cce_interlace:options','Finite positive integer reference_state and max_states required.');
end
if isfield(spec,'raw_only')&&(~islogical(spec.raw_only)||~isscalar(spec.raw_only))
    error('cce_interlace:options','raw_only must be logical.');
end
if isfield(spec,'condition_sites')&&~iscell(spec.condition_sites)
    error('cce_interlace:options','condition_sites must be a cell array.');
end
assignment=ensemble.assignments(spec.reference_state,:);
if numel(assignment)~=numel(model.bath_ids)||any(assignment<1)||...
   any(assignment>model.multiplicities(model.bath_ids))||any(assignment~=fix(assignment))
    error('cce_interlace:assignment','Interlacing requires a sampled product-state assignment.');
end
end

% Validate explicit choices and refuse unavailable extensions
function grumble(spin_system,parameters)
if ~isstruct(spin_system)||~isstruct(parameters)||~isscalar(parameters)
    error('spin_system and parameters must be scalar structures.');
end
required={'method','expansion','quantity','times','order',...
          'bath_state','meanfield','sequence','numerics'};
if ~all(isfield(parameters,required))||...
   (~isfield(parameters,'transition')&&~isfield(parameters,'rho0'))
    error('Required CCE parameter fields are missing.');
end
if ~ismember(parameters.method,{'cce','gcce','mecce','megcce','noise'})
    error('cce:capability','Method requires the cce_liouvillian or cce_noise extension.');
end
if ~ismember(parameters.expansion,{'standard','partition','disjoint','disjoint_inter'})
    error('cce:capability','Unsupported CCE expansion.');
end
if ~ismember(parameters.quantity,{'signal','density','observable','autocorrelation'})||...
   (ismember(parameters.method,{'cce','mecce'})&&~strcmp(parameters.quantity,'signal'))
    error('Projected CCE returns the selected-transition signal; gCCE also returns the observed-spin density matrix and observables.');
end
if strcmp(parameters.method,'noise')&&~strcmp(parameters.quantity,'autocorrelation')
    error('noise requires quantity=autocorrelation.');
end
if ~strcmp(parameters.meanfield,'none')&&~strcmp(parameters.meanfield,'static')
    error('meanfield must be explicitly none or static.');
end
if isfield(parameters,'averaging')&&...
   ~ismember(parameters.averaging,{'after_expansion','internal_then_expansion'})
    error('cce:capability','Interlaced averaging requires cce_interlace.');
end
if strcmp(parameters.bath_state.mode,'sampled')&&...
   ~isfield(parameters,'averaging')
    error('Sampled states require an explicit averaging estimator.');
end
for name={'drive','relaxation','convergence','control','pulse_family','time_axes'}
    if isfield(parameters,name{1})&&~isempty(parameters.(name{1}))
        error('cce:capability','Requested %s extension is not available in this baseline.',name{1});
    end
end
if isfield(parameters.numerics,'precision')&&~strcmp(parameters.numerics.precision,'double')
    error('cce:capability','Baseline propagation requires double precision.');
end
end


