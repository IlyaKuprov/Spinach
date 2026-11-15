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
% policy, physical ordering and elapsed time. This closed-system stage uses
% multiplicative reduction; density and additive noise are separate contracts.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)

% Check capabilities before constructing any retained Hilbert space
grumble(spin_system,pulse_sequence,parameters,assumptions);
started=tic; contract=cce_capabilities(pulse_sequence);
[parameters,policy]=cce_options(spin_system,parameters,pulse_sequence,assumptions);
if isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps)
    error('clexpn:dissipation','This experiment contract is closed-system only.');
end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    error('clexpn:density','Density requires the additive density contract.');
end
if ~strcmp(parameters.cce.expansion,'standard')||...
   ~strcmp(parameters.cce.bath_state.mode,'mixed_product')&&...
   ~strcmp(parameters.cce.bath_state.mode,'explicit_product')
    error('clexpn:ensemble','This stage supports standard expansion of explicit or mixed product states.');
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

% Reject policy controls whose context implementation belongs to the next stage
if isfield(parameters.cce,'selection')&&isfield(parameters.cce.selection,'adaptive')
    error('clexpn:adaptive','Adaptive convergence policy is not yet available.');
end

% Compile immutable physical data, not a full-bath basis
model=cce_compile(spin_system,parameters,assumptions);
if ~ismember(model.formalism,contract.formalisms)
    error('clexpn:formalism','Experiment does not support the requested formalism.');
end
if isempty(parameters.cce.order)
    parameters.cce.order=numel(model.bath_ids);
    policy.defaults{end+1}='parameters.cce.order';
end

% Existing graph and cluster services receive their own algorithm options
options=parameters.cce;
selection=struct('graph',options.graph,'order',options.order,'numerics',options.numerics);
selection.method='gcce'; selection.expansion=options.expansion;
if strcmp(options.approximation,'projected'), selection.method='cce'; end
for field={'groups','selection'}
    if isfield(options,field{1}), selection.(field{1})=options.(field{1}); end
end
graph=cce_graph(model,selection); clusters=cce_clusters(model,graph,selection);
configuration=struct('configuration_id',1,'state_id',1,'weight',1);
ensemble.spin_config=configuration; ensemble.weights=1;
schedule.pulse_sequence=pulse_sequence; schedule.contract=contract;

% Scheduler calls the registered evaluator, which calls the five-argument experiment
raw=cce_schedule(model,clusters,ensemble,schedule,parameters,contract.evaluator);
reduction=struct('method','gcce','numerics',options.numerics,'reducer',options.reducer);
reduced=cce_reduce(raw{1},clusters,reduction); answer=reduced.value;

% Keep diagnostics separate from the natural sequence output
diagnostics.raw=raw; diagnostics.clusters=clusters; diagnostics.graph=graph;
diagnostics.reduction=reduced; diagnostics.policy=policy;
diagnostics.parameters=parameters; diagnostics.contract=contract;
diagnostics.physical_ids=model.physical_ids; diagnostics.observed_ids=model.observed_ids;
diagnostics.bath_ids=model.bath_ids; diagnostics.elapsed_seconds=toc(started);

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


