% Cluster correlation expansion for a prepared Spinach spin system. Syntax:
%
%                       result=cce(spin_system,parameters)
%
% Required fields: central_spins, bath_spins, method (cce/gcce/mecce/megcce/noise), expansion
% (standard), quantity (coherence/density/observable/autocorrelation), times (whole-sequence
% seconds), order, central, assumptions, orientation, interaction_cutoff,
% graph, bath_state, meanfield, sequence, numerics.atol/rtol/max_memory_bytes.
% cce_prepare, cce_graph, cce_states and cce_sequence document their schemas.
% normalized=true divides gcce observables by the supplied initial state;
% conditional branch coherence is intrinsically normalized. Density is not.
% reducer is powers (default) or direct. Sampled estimators explicitly set
% averaging='after_expansion' or 'internal_then_expansion'.
%
% Outputs include complex value/coherence, per-order estimates, validity,
% reasons, propagated error bounds, raw clusters and exact-full-bath bypass
% metadata. Density outputs additionally report physicality without repair.
% Extensions fail with cce:capability until their implementations are present;
% no optional control, checkpoint or dissipative physics is silently ignored.
% ME jumps use ordered physical support, dimensionless operator and s^-1 rate
% fields (cce_liouvillian). Noise requires explicit ordering, connected flag
% and Hermitian support/operator terms (cce_noise); its reducer is additive.
%
% ilya.kuprov@weizmann.ac.il

function result=cce(spin_system,parameters)

% Validate orchestration choices before canonicalizing the physical model
grumble(spin_system,parameters);
started=tic; model=cce_prepare(spin_system,parameters);

% Dispatch open-system and additive-noise cluster evaluators
evaluator=@cce_cluster;
if ismember(parameters.method,{'mecce','megcce'}), evaluator=@cce_master; end
if strcmp(parameters.method,'noise'), evaluator=@cce_noise; end
state_parameters=parameters;
if ismember(parameters.method,{'mecce','noise'}), state_parameters.method='cce'; end
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

% Apply W12 conditional averages before irreducible subset reduction
if strcmp(parameters.bath_state.mode,'interlaced')
    if ~strcmp(parameters.quantity,'coherence')
        error('cce:interlaced','Interlaced multiplication requires scalar coherence.');
    end
    empty=struct('id',uint64(0),'sites',[]);
    estimates=zeros(numel(parameters.times),numel(ensemble.contexts));
    by_order=zeros(numel(parameters.times),parameters.order,numel(ensemble.contexts));
    valid=true(size(estimates)); bounds=zeros(size(estimates)); records=cell(numel(ensemble.contexts),numel(clusters));
    for s=1:numel(ensemble.contexts)
        conditioned=parameters; conditioned.interlaced.reference_state=s;
        baseline=cce_interlace(model,empty,ensemble,schedule,conditioned);
        product=baseline.value; uncertainty=baseline.abs_error_bound;
        for order=1:parameters.order
            for j=find([clusters.order]==order)
                item=cce_interlace(model,clusters(j),ensemble,schedule,conditioned);
                uncertainty=abs(product).*item.factor_abs_error_bound+...
                    abs(item.factor).*uncertainty+uncertainty.*item.factor_abs_error_bound;
                product=product.*item.factor; valid(:,s)=valid(:,s)&item.factor_valid;
                records{s,j}=item;
            end
            by_order(:,order,s)=product;
        end
        estimates(:,s)=product; bounds(:,s)=uncertainty;
    end
    value=estimates*ensemble.weights(:); uncertainty=bounds*ensemble.weights(:);
    mask=all(valid,2)&uncertainty<=parameters.numerics.atol+parameters.numerics.rtol*abs(value);
    result=struct('value',value,'coherence',value,'valid',mask,'times',parameters.times,...
                  'abs_error_bound',uncertainty,...
                  'by_order',sum(by_order.*reshape(ensemble.weights,1,1,[]),3),...
                  'raw',{records},'clusters',clusters,'graph',graph,'schedule',schedule,...
                  'ensemble',ensemble,'model_metadata',rmfield(model,'spin_system'),...
                  'status','complete','diagnostics',struct('conditioning','W12_shared_K'));
    result.reason=repmat({''},size(value)); result.reason(~mask)={'unresolved_interlaced_factor'};
    result.elapsed=toc(started); result.sampling_standard_error=NaN(size(value));
    if numel(ensemble.contexts)>1
        result.sampling_standard_error=std(estimates,0,2)/sqrt(numel(ensemble.contexts));
    end
    if any(~mask), result.status='incomplete'; end
    return;
end

% Propagate each globally assigned bath state before ensemble reduction
ns=numel(ensemble.contexts); reductions=cell(1,ns);
all_raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator);
full_index=find(arrayfun(@(c) isequal(sort(c.sites),sort(model.bath_ids)),clusters),1);
if isempty(model.bath_ids), full_index=0; end

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
if ns>1&&strcmp(averaging,'after_expansion')
    for s=1:ns
        result.sampling_standard_error=result.sampling_standard_error+...
            abs(reductions{s}.value-result.value).^2/(ns*(ns-1));
    end
    result.sampling_standard_error=sqrt(result.sampling_standard_error);
end
result.coherence=result.value; result.raw=all_raw; result.times=parameters.times;
result.clusters=clusters; result.graph=graph; result.schedule=schedule;
result.partition=partition;
result.ensemble=ensemble; result.model_hash=model.source_hash;
result.input_parameters=parameters; result.averaging=averaging;
result.method=parameters.method; result.exact_reference=~isempty(full_index);
if result.exact_reference, result.method='exact-reference'; end
result.elapsed_seconds=toc(started); result.convergence='not_requested';
result.budget_satisfied=result.valid&...
    result.abs_error_bound<=parameters.numerics.atol+parameters.numerics.rtol*abs(result.value);

% Diagnose a reconstructed central density without symmetrization or clipping
if strcmp(parameters.quantity,'density')
    d=ensemble.central_dimension; nt=numel(parameters.times);
    result.density=reshape(result.value.',d,d,nt);
end
result.diagnostics=cce_diagnostics(model,all_raw,result,parameters);

end

% Validate explicit choices and refuse unavailable extensions
function grumble(spin_system,parameters)
if ~isstruct(spin_system)||~isstruct(parameters)||~isscalar(parameters)
    error('spin_system and parameters must be scalar structures.');
end
required={'method','expansion','quantity','times','order','central',...
          'bath_state','meanfield','sequence','numerics'};
if ~all(isfield(parameters,required)), error('Required CCE parameter fields are missing.'); end
if ~ismember(parameters.method,{'cce','gcce','mecce','megcce','noise'})
    error('cce:capability','Method requires the cce_liouvillian or cce_noise extension.');
end
if ~ismember(parameters.expansion,{'standard','partition','disjoint','disjoint_inter'})
    error('cce:capability','Unsupported CCE expansion.');
end
if ~ismember(parameters.quantity,{'coherence','density','observable','autocorrelation'})||...
   (ismember(parameters.method,{'cce','mecce'})&&~strcmp(parameters.quantity,'coherence'))
    error('Conditional CCE returns branch coherence; gcce also returns density and observables.');
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

