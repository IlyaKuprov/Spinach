% Cluster-expansion context for registered magnetic resonance experiments.
% Syntax:
%
% [answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)
%
% spin_system is create output; no full-bath basis is required. pulse_sequence
% is a registered experiment with signature (spin_system,parameters,H,R,K).
% parameters.observed contains physical spin IDs; times contains elapsed
% sequence times in seconds. transition or explicit rho0/coil defines the
% observed transition or observable. Algorithm controls belong to cce.
% assumptions is the Spinach interaction assumption string.
%
% answer is the experiment's natural numeric output, without a legacy result
% wrapper. diagnostics includes raw cluster responses, reduction, resolved
% policy, physical ordering, and elapsed time. This closed-system stage uses
% multiplicative reduction; density and additive noise are separate contracts.
% For noise_correlation, parameters.spectrum requests cce_spectrum after the
% completed ensemble reduction; the result is diagnostics.spectrum.
%
% Parameters:
%
%   spin_system    - native Spinach system returned by create
%   pulse_sequence - registered five-argument experiment function handle
%   parameters     - experiment parameters and nested cce algorithm controls
%   assumptions    - explicit Spinach Hamiltonian assumption string
%
% Supply parameters.orientation as a 1-by-3 Euler row in radians (default
% [0 0 0]), or parameters.grid as a grid name or explicit alphas, betas,
% gammas, and weights vectors. Positive weights are normalised; zero-weight
% nodes are omitted in both output modes. Each node has its own physical
% model and cluster generators; the receiver and spin-state draws stay fixed.
% A checkpoint requires exactly one retained node. Grid nodes run serially;
% cluster/state parallelism stays with the scheduler, or with an outer worker.
%
% On a grid, parameters.sum_up=true (default) averages completed expansions
% unless cce.averaging='internal_then_expansion', which averages raw cluster
% records and observed-only references over nodes and states before reduction.
% Internal averaging requires a common physical cluster family. sum_up=false
% returns unweighted completed outputs as a column cell in input grid order.
% Without a grid, the single-orientation answer is always numeric.
%
% Outputs:
%
%   answer      - experiment output on its declared time/observable axes
%   diagnostics - physical ordering, cluster records, reduction, resolved
%                 policy, ensemble, and elapsed execution time. Grid output
%                 has grid (evaluated nodes and normalised weights), grid_indices
%                 (original one-based indices), and matching column cells values
%                 and orientations. Only a summed grid has a combined reduction.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)

% Check capabilities before constructing any retained Hilbert space
grumble(spin_system,pulse_sequence,parameters,assumptions);
started=tic;
contract=cce_capabilities(pulse_sequence);
[parameters,policy]=cce_options(spin_system,parameters,pulse_sequence,assumptions);
if strcmp(parameters.cce.approximation,'projected')&&~contract.projected
    error('clexpn:projection','This experiment requires generalized evolution.');
end
if (isfield(parameters,'noise')||isfield(parameters,'spectrum'))&&...
   ~isequal(pulse_sequence,@noise_correlation)
    error('clexpn:noise','Noise and spectrum parameters belong to the noise_correlation experiment.');
end
if isequal(pulse_sequence,@noise_correlation)&&...
   strcmp(parameters.cce.bath_state.mode,'interlaced')
    error('clexpn:noise','Connected noise uses ordinary ensemble moment averaging, not interlaced factors.');
end
if isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps)&&~contract.dissipative
    error('clexpn:dissipation','This experiment contract is closed-system only.');
end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')&&~contract.density_reduction
    error('clexpn:density','Density requires a registered density output contract.');
end

% Derive default numerical controls from precision and available host memory
numerics=parameters.cce.numerics;
if ~isfield(numerics,'atol')
    numerics.atol=sqrt(eps);
    policy.defaults{end+1}='parameters.cce.numerics.atol';
end
if ~isfield(numerics,'rtol')
    numerics.rtol=sqrt(eps);
    policy.defaults{end+1}='parameters.cce.numerics.rtol';
end
if ~isfield(numerics,'max_memory_bytes')
    operating_system=java.lang.management.ManagementFactory.getOperatingSystemMXBean();
    numerics.max_memory_bytes=floor(double(operating_system.getFreeMemorySize())/2);
    policy.defaults{end+1}='parameters.cce.numerics.max_memory_bytes';
end
parameters.cce.numerics=numerics;

% Reject adaptive controllers until the convergence stage owns their execution
if isfield(parameters.cce,'selection')&&isfield(parameters.cce.selection,'adaptive')
    error('clexpn:adaptive','Adaptive convergence policy is not yet available.');
end

% Resolve the quadrature before any checkpoint or cluster evaluation
if isfield(parameters,'grid')
    grid=parameters.grid;
    indices=find(grid.weights(:)>0);
    if numel(indices)>1&&isfield(parameters.cce,'checkpoint')
        error('clexpn:checkpoint','Checkpoint files require exactly one positive-weight orientation.');
    end
    for field={'alphas','betas','gammas','weights'}
        data=grid.(field{1})(:);
        grid.(field{1})=double(data(indices));
    end
    grid.weights=grid.weights/sum(grid.weights);
    if ~isempty(getCurrentTask())
        parameters.cce.parallel.mode='serial';
    end
end

% Use the supplied finite bath when no expansion order was declared
if isempty(parameters.cce.order)
    parameters.cce.order=numel(parameters.cce.bath_spins);
    policy.defaults{end+1}='parameters.cce.order';
end

% Disjoint families retain single groups or all inter-group pairs
if strcmp(parameters.cce.expansion,'disjoint')
    parameters.cce.order=1;
end
if strcmp(parameters.cce.expansion,'disjoint_inter')
    parameters.cce.order=2;
end

% Construct product-state assignments without changing the caller's random stream
state=parameters.cce.bath_state;
ids=parameters.cce.bath_spins;
count=1;
if ismember(state.mode,{'sampled','interlaced'})
    count=state.samples;
end
bytes=count*(1024+16*sum(spin_system.comp.mults(ids).^2)+32*numel(ids));
if ~isfinite(bytes)||bytes>numerics.max_memory_bytes
    error('clexpn:memory','Sampled product-state storage exceeds the declared memory budget.');
end
configuration=struct('configuration_id',1,'state_id',1,'weight',1/count,...
                     'bath_ids',ids,'bath_rho',{{}},'assignment',zeros(1,numel(ids)));
if isfield(state,'configuration_id')
    configuration.configuration_id=state.configuration_id;
end
configs=repmat(configuration,1,count);
bases=cell(1,numel(ids));
probabilities=bases;
for n=1:count
    if ismember(state.mode,{'sampled','interlaced'})
        stream=RandStream('Threefry','Seed',state.seed);
        stream.Substream=(configuration.configuration_id-1)*count+n;
    end
    configs(n).state_id=n;
    configs(n).bath_rho=cell(1,numel(ids));
    for j=1:numel(ids)
        d=spin_system.comp.mults(ids(j));
        switch state.mode
            case 'mixed_product'
                rho=eye(d)/d;
            case 'explicit_product'
                rho=state.rho{j};
            otherwise
                if ischar(state.basis)
                    bases{j}=eye(d);
                else
                    bases{j}=state.basis{j};
                end
                if ischar(state.probabilities)
                    probabilities{j}=ones(d,1)/d;
                else
                    probabilities{j}=state.probabilities{j}(:);
                end
                u=rand(stream);
                pick=find(u<=cumsum(probabilities{j}),1);
                if isempty(pick)
                    pick=d;
                end
                configs(n).assignment(j)=pick;
                vector=bases{j}(:,pick);
                rho=vector*vector';
        end
        configs(n).bath_rho{j}=rho;
    end
end
ensemble=struct('spin_config',configs,'weights',[configs.weight],...
                'bases',{bases},'probabilities',{probabilities});
schedule.pulse_sequence=pulse_sequence;
schedule.contract=contract;

% A single crystal returns the experiment's natural numeric output
if ~isfield(parameters,'grid')
    if isfield(parameters,'sum_up')
        parameters=rmfield(parameters,'sum_up');
    end
    diagnostics=crystallite(spin_system,schedule,parameters,assumptions,policy,ensemble);
    answer=diagnostics.reduction.value;
    diagnostics.elapsed_seconds=toc(started);
    return;
end

% Rotate the physical model independently at each positive-weight grid node
local=rmfield(parameters,{'grid','sum_up'});
count=numel(indices);
orientations=cell(count,1);
values=cell(count,1);
for n=1:count
    local.orientation=[grid.alphas(n) grid.betas(n) grid.gammas(n)];
    orientations{n}=crystallite(spin_system,schedule,local,assumptions,policy,ensemble);
    values{n}=orientations{n}.reduction.value;
end

% Expose the original grid indices alongside the unweighted, ordered results
diagnostics.grid=grid;
diagnostics.grid_indices=indices;
diagnostics.values=values;
diagnostics.orientations=orientations;
diagnostics.parameters=parameters;
if ~parameters.sum_up
    answer=values;
    diagnostics.elapsed_seconds=toc(started);
    return;
end

% Choose where quadrature enters the cluster expansion
first=orientations{1};
reduction=first.reduction;
if count>1&&strcmp(parameters.cce.averaging,'internal_then_expansion')

    % Pool raw state records, including their observed-only references
    samples=numel(first.raw);
    raw=cell(1,count*samples);
    for n=1:count
        current=orientations{n};
        if ~isequal({current.clusters.sites},{first.clusters.sites})||...
           ~isequal({current.clusters.groups},{first.clusters.groups})||...
           ~isequal(current.graph.bath_ids,first.graph.bath_ids)
            error('clexpn:clusters','Internal averaging requires identical physical cluster families.');
        end
        for k=1:samples
            records=current.raw{k};
            for j=1:numel(records)
                records(j).weight=grid.weights(n)*records(j).weight;
            end
            raw{(n-1)*samples+k}=records;
        end
    end
    request=first.parameters;
    request.cce.bath_spins=first.graph.bath_ids;
    reduction=cce_reduce(raw,first.clusters,request);

    % Only complete, uncentred records have an internal-first sampling estimate
    reduction.sampling_standard_error=first.reduction.sampling_standard_error;
    complete=reduction.complete_record;
    if strcmp(parameters.cce.bath_state.mode,'sampled')&&samples>1&&...
       ~isempty(complete)&&~(isfield(parameters,'noise')&&parameters.noise.connected)
        draws=zeros([size(reduction.value) samples]);
        for n=1:count
            for k=1:samples
                draws(:,:,k)=draws(:,:,k)+grid.weights(n)*orientations{n}.raw{k}(complete).value;
            end
        end
        reduction.sampling_standard_error=std(draws,0,3)/sqrt(samples);
    end

elseif count>1

    % Combine completed expansions and error ledgers, then each shared spin draw
    parts=cellfun(@(node)node.reduction,orientations,'UniformOutput',false);
    reduction=weighted_reduction(parts,grid.weights);
    samples=numel(reduction.samples);
    for k=1:samples
        draws=cellfun(@(part)part.samples{k},parts,'UniformOutput',false);
        reduction.samples{k}=weighted_reduction(draws,grid.weights);
    end
    if samples>1
        draws=cellfun(@(part)part.value,reduction.samples,'UniformOutput',false);
        reduction.sampling_standard_error=std(cat(3,draws{:}),0,3)/sqrt(samples);
    end

end
answer=reduction.value;
diagnostics.reduction=reduction;
if isfield(parameters,'spectrum')
    diagnostics.spectrum=cce_spectrum(parameters.times,answer,struct('spectrum',parameters.spectrum));
end
diagnostics.elapsed_seconds=toc(started);

end

% Evaluate the experiment using only this crystallite's retained cluster generators
function diagnostics=crystallite(spin_system,schedule,parameters,assumptions,policy,ensemble)

% Use the registered experiment and its validated output contract
started=tic;
pulse_sequence=schedule.pulse_sequence;
contract=schedule.contract;

% Compile immutable physical data, not a full-bath basis
model=cce_compile(spin_system,parameters,assumptions);
model.pulse_sequence=pulse_sequence;
if isfield(model,'transition')
    policy.transition=model.transition;
    policy.transition_meta=model.transition_meta;
end
if isfield(parameters.cce,'basis')&&isfield(parameters.cce.basis,'formalism')
    model.formalism=parameters.cce.basis.formalism;
elseif isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps)
    model.formalism='zeeman-liouv';
    policy.defaults{end+1}='parameters.cce.basis.formalism';
end
if ~ismember(model.formalism,contract.formalisms)
    error('clexpn:formalism','The experiment does not support the requested representation.');
end

% Validate propagation choices here rather than in graph or reduction services
if isfield(parameters.cce.numerics,'propagator')
    supported={'eig'};
    if strcmp(model.formalism,'zeeman-liouv')
        supported={'stock'};
    elseif isequal(pulse_sequence,@hahn_echo)&&strcmp(parameters.cce.approximation,'generalized')
        supported={'eig','expm'};
    end
    if ~ischar(parameters.cce.numerics.propagator)||~isrow(parameters.cce.numerics.propagator)||...
       ~ismember(parameters.cce.numerics.propagator,supported)
        error('clexpn:propagator','Propagation preference is incompatible with this experiment and representation.');
    end
end

% Build the physical connectivity and downward-closed cluster family
options=parameters.cce;
selection=struct('graph',options.graph,'order',options.order,'numerics',options.numerics);
if isfield(selection.numerics,'propagator')
    selection.numerics=rmfield(selection.numerics,'propagator');
end
selection.method='gcce';
selection.expansion=options.expansion;
if strcmp(options.approximation,'projected')
    selection.method='cce';
end
if isequal(pulse_sequence,@noise_correlation)
    selection.method='noise';
    selection.noise=parameters.noise;
    selection.bath_state=options.bath_state;
end
if strcmp(options.approximation,'projected')
    selection.second_order=options.second_order;
end
for field={'groups','selection','bath_radius'}
    if isfield(options,field{1})
        selection.(field{1})=options.(field{1});
    end
end
graph=cce_graph(model,selection);
if isfield(options,'partition')
    partition=cce_partition(model,graph,parameters);
    if ~partition.complete
        error('clexpn:partition','Partition adaptation is incomplete: %s.',partition.diagnostics.reason);
    end
    selection.groups=partition.groups;
    if ~isempty(partition.groups)
        parameters.cce.groups=partition.groups;
    end
    model.parameters=parameters;
    model.groups=partition.groups;
    policy.partition=partition;
end
clusters=cce_clusters(model,graph,selection);
model.parameters=parameters;

% Scheduler calls the registered evaluator, which calls the five-argument experiment
raw=cce_schedule(model,clusters,ensemble,schedule,parameters,contract.evaluator);

% Test completeness on the graph domain without trimming the exterior state
request=parameters;
request.cce.bath_spins=graph.bath_ids;
reduced=cce_reduce(raw,clusters,request);

% Keep diagnostics separate from the natural sequence output
diagnostics.raw=raw;
diagnostics.clusters=clusters;
diagnostics.graph=graph;
diagnostics.reduction=reduced;
diagnostics.policy=policy;
diagnostics.parameters=parameters;
diagnostics.contract=contract;
diagnostics.physical_ids=model.physical_ids;
diagnostics.observed_ids=model.observed_ids;
diagnostics.bath_ids=model.bath_ids;
diagnostics.ensemble=ensemble;
if isfield(parameters,'spectrum')
    diagnostics.spectrum=cce_spectrum(parameters.times,reduced.value,...
                                     struct('spectrum',parameters.spectrum));
end
diagnostics.elapsed_seconds=toc(started);

end

% Integrate completed values and absolute bounds without concealing failed nodes
function reduced=weighted_reduction(parts,weights)

% Each validity mask and reason accompanies the same output element and order
reduced=parts{1};
for field={'value','by_order','abs_error_bound','by_order_error'}
    name=field{1};
    reduced.(name)=weights(1)*parts{1}.(name);
    for n=2:numel(parts)
        reduced.(name)=reduced.(name)+weights(n)*parts{n}.(name);
    end
end
for n=2:numel(parts)
    current=parts{n};
    reduced.valid=reduced.valid&current.valid;
    reduced.reason(~current.valid)=current.reason(~current.valid);
    reduced.by_order_valid=reduced.by_order_valid&current.by_order_valid;
    failed=~current.by_order_valid;
    reduced.by_order_reason(failed)=current.by_order_reason(failed);
end

end

% Validate the public context signature before any expensive work
function grumble(spin_system,pulse_sequence,parameters,assumptions)
if ~isstruct(spin_system)||~isfield(spin_system,'comp')
    error('clexpn:spin_system','spin_system must be create output.');
end
if ~isa(pulse_sequence,'function_handle')||~isstruct(parameters)||...
   ~isscalar(parameters)||~ischar(assumptions)||~isrow(assumptions)
    error('clexpn:input','Supply an experiment handle, parameter structure and assumptions string.');
end
end


