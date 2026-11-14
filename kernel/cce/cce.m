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
graph=cce_graph(model,parameters); clusters=cce_clusters(model,graph,parameters);

% Propagate each globally assigned bath state before ensemble reduction
ns=numel(ensemble.contexts); all_raw=cell(1,ns); reductions=cell(1,ns);
empty=struct('id',uint64(0),'sites',[]);
full_index=find(arrayfun(@(c) isequal(sort(c.sites),sort(model.bath_ids)),clusters),1);
if isempty(model.bath_ids), full_index=0; end
for s=1:ns
    records=evaluator(model,empty,ensemble.contexts(s),schedule,parameters);
    for k=1:numel(clusters)
        records(k+1)=evaluator(model,clusters(k),ensemble.contexts(s),schedule,parameters);
    end
    all_raw{s}=records;
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
if ns>1&&strcmp(averaging,'after_expansion')
    for s=1:ns
        result.sampling_standard_error=result.sampling_standard_error+...
            abs(reductions{s}.value-result.value).^2/(ns*(ns-1));
    end
    result.sampling_standard_error=sqrt(result.sampling_standard_error);
end
result.coherence=result.value; result.raw=all_raw; result.times=parameters.times;
result.clusters=clusters; result.graph=graph; result.schedule=schedule;
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
    result.diagnostics.trace_error=nan(nt,1);
    result.diagnostics.hermiticity_error=nan(nt,1);
    result.diagnostics.minimum_eigenvalue=nan(nt,1);
    for k=1:nt
        rho=result.density(:,:,k);
        if all(isfinite(rho),'all')
            result.diagnostics.trace_error(k)=abs(trace(rho)-1);
            result.diagnostics.hermiticity_error(k)=norm(rho-rho','fro');
            result.diagnostics.minimum_eigenvalue(k)=min(eig((rho+rho')/2));
        end
    end
end

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
if ~strcmp(parameters.expansion,'standard')
    error('cce:capability','Nonstandard grouping requires the cce_partition extension.');
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
if ismember(parameters.bath_state.mode,{'sampled','interlaced'})&&...
   ~isfield(parameters,'averaging')
    error('Sampled states require an explicit averaging estimator.');
end
if strcmp(parameters.bath_state.mode,'interlaced')
    error('cce:capability','Interlaced states require cce_interlace.');
end
for name={'drive','relaxation','convergence','checkpoint','control','pulse_family','time_axes'}
    if isfield(parameters,name{1})&&~isempty(parameters.(name{1}))
        error('cce:capability','Requested %s extension is not available in this baseline.',name{1});
    end
end
if isfield(parameters,'parallel')&&~isempty(parameters.parallel)
    if ~isstruct(parameters.parallel)||~isfield(parameters.parallel,'mode')||...
       ~strcmp(parameters.parallel.mode,'serial')
        error('cce:capability','Nonserial execution requires cce_schedule.');
    end
end
if isfield(parameters.numerics,'precision')&&~strcmp(parameters.numerics.precision,'double')
    error('cce:capability','Baseline propagation requires double precision.');
end
if isfield(parameters,'deer')&&isfield(parameters.deer,'enabled')&&parameters.deer.enabled
    error('cce:capability','DEER selection requires cce_controls.');
end
end

