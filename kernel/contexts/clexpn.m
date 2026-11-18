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
% Outputs:
%
%   answer      - experiment output on its declared time/observable axes
%   diagnostics - physical ordering, cluster records, reduction, resolved
%                 policy, ensemble, and elapsed execution time
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)

% Check capabilities before constructing any retained Hilbert space
grumble(spin_system,pulse_sequence,parameters,assumptions);
started=tic; contract=cce_capabilities(pulse_sequence);
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
    numerics.atol=sqrt(eps); policy.defaults{end+1}='parameters.cce.numerics.atol';
end
if ~isfield(numerics,'rtol')
    numerics.rtol=sqrt(eps); policy.defaults{end+1}='parameters.cce.numerics.rtol';
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

% Compile immutable physical data, not a full-bath basis
model=cce_compile(spin_system,parameters,assumptions);
model.pulse_sequence=pulse_sequence;
if isfield(model,'transition')
    policy.transition=model.transition; policy.transition_meta=model.transition_meta;
end
if isfield(parameters.cce,'basis')&&isfield(parameters.cce.basis,'formalism')
    model.formalism=parameters.cce.basis.formalism;
elseif isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps)
    model.formalism='zeeman-liouv'; policy.defaults{end+1}='parameters.cce.basis.formalism';
end
if ~ismember(model.formalism,contract.formalisms)
    error('clexpn:formalism','The experiment does not support the requested representation.');
end

% Validate propagation choices here rather than in graph or reduction services
if isfield(numerics,'propagator')
    supported={'eig'};
    if strcmp(model.formalism,'zeeman-liouv')
        supported={'stock'};
    elseif isequal(pulse_sequence,@hahn_echo)&&strcmp(parameters.cce.approximation,'generalized')
        supported={'eig','expm'};
    end
    if ~ischar(numerics.propagator)||~isrow(numerics.propagator)||...
       ~ismember(numerics.propagator,supported)
        error('clexpn:propagator','Propagation preference is incompatible with this experiment and representation.');
    end
end

% Use the supplied finite bath when no expansion order was declared
if isempty(parameters.cce.order)
    parameters.cce.order=numel(model.bath_ids);
    policy.defaults{end+1}='parameters.cce.order';
end

% Disjoint families retain single groups or all inter-group pairs
if strcmp(parameters.cce.expansion,'disjoint'), parameters.cce.order=1; end
if strcmp(parameters.cce.expansion,'disjoint_inter'), parameters.cce.order=2; end

% Build the physical connectivity and downward-closed cluster family
options=parameters.cce;
selection=struct('graph',options.graph,'order',options.order,'numerics',options.numerics);
if isfield(selection.numerics,'propagator')
    selection.numerics=rmfield(selection.numerics,'propagator');
end
selection.method='gcce'; selection.expansion=options.expansion;
if strcmp(options.approximation,'projected'), selection.method='cce'; end
if isequal(pulse_sequence,@noise_correlation)
    selection.method='noise'; selection.noise=parameters.noise;
    selection.bath_state=options.bath_state;
end
if strcmp(options.approximation,'projected'), selection.second_order=options.second_order; end
for field={'groups','selection','bath_radius'}
    if isfield(options,field{1}), selection.(field{1})=options.(field{1}); end
end
graph=cce_graph(model,selection);
if isfield(options,'partition')
    partition=cce_partition(model,graph,parameters);
    if ~partition.complete
        error('clexpn:partition','Partition adaptation is incomplete: %s.',partition.diagnostics.reason);
    end
    selection.groups=partition.groups;
    if ~isempty(partition.groups), parameters.cce.groups=partition.groups; end
    model.parameters=parameters; model.groups=partition.groups; policy.partition=partition;
end
clusters=cce_clusters(model,graph,selection);
model.parameters=parameters;

% Construct product-state assignments without changing the caller's random stream
state=parameters.cce.bath_state;
ids=model.bath_ids; count=1;
if ismember(state.mode,{'sampled','interlaced'}), count=state.samples; end
bytes=count*(1024+16*sum(model.multiplicities(ids).^2)+32*numel(ids));
if ~isfinite(bytes)||bytes>model.numerics.max_memory_bytes
    error('clexpn:memory','Sampled product-state storage exceeds the declared memory budget.');
end
configuration=struct('configuration_id',1,'state_id',1,'weight',1/count,...
                     'bath_ids',ids,'bath_rho',{{}},'assignment',zeros(1,numel(ids)));
if isfield(state,'configuration_id'), configuration.configuration_id=state.configuration_id; end
configs=repmat(configuration,1,count); bases=cell(1,numel(ids)); probabilities=bases;
for n=1:count
    if ismember(state.mode,{'sampled','interlaced'})
        stream=RandStream('Threefry','Seed',state.seed);
        stream.Substream=(configuration.configuration_id-1)*count+n;
    end
    configs(n).state_id=n; configs(n).bath_rho=cell(1,numel(ids));
    for j=1:numel(ids)
        d=model.multiplicities(ids(j));
        switch state.mode
            case 'mixed_product'
                rho=eye(d)/d;
            case 'explicit_product'
                rho=state.rho{j};
            otherwise
                if ischar(state.basis), bases{j}=eye(d); else, bases{j}=state.basis{j}; end
                if ischar(state.probabilities)
                    probabilities{j}=ones(d,1)/d;
                else
                    probabilities{j}=state.probabilities{j}(:);
                end
                u=rand(stream); pick=find(u<=cumsum(probabilities{j}),1);
                if isempty(pick), pick=d; end
                configs(n).assignment(j)=pick;
                vector=bases{j}(:,pick); rho=vector*vector';
        end
        configs(n).bath_rho{j}=rho;
    end
end
ensemble=struct('spin_config',configs,'weights',[configs.weight],...
                'bases',{bases},'probabilities',{probabilities});
schedule.pulse_sequence=pulse_sequence; schedule.contract=contract;

% Scheduler calls the registered evaluator, which calls the five-argument experiment
raw=cce_schedule(model,clusters,ensemble,schedule,parameters,contract.evaluator);
reduced=cce_reduce(raw,clusters,parameters); answer=reduced.value;

% Keep diagnostics separate from the natural sequence output
diagnostics.raw=raw; diagnostics.clusters=clusters; diagnostics.graph=graph;
diagnostics.reduction=reduced; diagnostics.policy=policy;
diagnostics.parameters=parameters; diagnostics.contract=contract;
diagnostics.physical_ids=model.physical_ids; diagnostics.observed_ids=model.observed_ids;
diagnostics.bath_ids=model.bath_ids; diagnostics.ensemble=ensemble;
if isfield(parameters,'spectrum')
    diagnostics.spectrum=cce_spectrum(parameters.times,answer,...
                                     struct('spectrum',parameters.spectrum));
end
diagnostics.elapsed_seconds=toc(started);

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


