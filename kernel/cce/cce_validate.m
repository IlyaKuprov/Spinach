% Reject unsupported CCE parameter fields before any evaluator dispatch.
% Syntax: cce_validate(parameters)
%
% This shared contract checks optional field names and method applicability;
% individual helpers validate the dimensions and values they consume.
%
% ilya.kuprov@weizmann.ac.il

function cce_validate(parameters)

% Validate the public container before inspecting optional fields
grumble(parameters);

% Reject misspelled top-level parameters consistently in standalone helpers
known={'central_spins','bath_spins','method','expansion','quantity','times',...
       'order','central','assumptions','orientation','hamiltonian_backend',...
       'interaction_cutoff','bath_radius','graph','groups','partition',...
       'bath_state','sampling','meanfield','meanfield_radius','sequence',...
       'controls','deer','jumps','noise','spectrum','numerics','averaging',...
       'analysis','diagnostics','run','checkpoint','retain_raw','normalized',...
       'operators','interlaced','convergence','reducer','parallel','initial_rho',...
       'clusters','cluster_cap','second_order'};
fields(parameters,known,'parameters');

% Reject the obsolete control container even when it is empty
if isfield(parameters,'controls')
    error('cce:controls','parameters.controls is unsupported; use sequence events and parameters.deer.');
end

% Reject reserved driver extensions rather than silently dropping them
for name={'sampling','run','diagnostics'}
    if isfield(parameters,name{1})
        error('cce:capability','parameters.%s is unsupported; use the corresponding helper explicitly.',name{1});
    end
end

% Check the implemented nested parameter vocabulary
contracts={'central',{'alpha','beta','rho'};...
           'second_order',{'enabled','gap_tolerance','overlap_tolerance'};...
           'numerics',{'atol','rtol','max_memory_bytes','precision','zero_policy','propagator'};...
           'graph',{'mode','score','radius','threshold','adjacency','scores','mismatch','mismatch_scale','approximate'};...
           'bath_state',{'mode','rho','basis','probabilities','samples','seed','configuration_id'};...
           'parallel',{'mode','batch_size','max_tasks'};...
           'deer',{'enabled','detuning_center','halfwidth','species','sites','allow_empty'};...
           'checkpoint',{'file','resume'};...
           'noise',{'terms','connected','ordering','branch'};...
           'interlaced',{'condition_sites','reference_state','max_states','raw_only'};...
           'partition',{'mode','max_dimension','group_count','adaptive','patches','patch_atol','patch_rtol','pilot_contexts','heldout_contexts','max_rounds','max_probes'};...
           'spectrum',{'frequency_grid','convention','window','weights','sidedness','time_domain','stationary','modulation'};...
           'analysis',{'observable','envelope','fit','model','interval','weights','confidence','confidence_level'}};
for k=1:size(contracts,1)
    name=contracts{k,1};
    if isfield(parameters,name), fields(parameters.(name),contracts{k,2},['parameters.' name]); end
end

% Distinguish named controls from explicit instantaneous and finite events
if isfield(parameters,'sequence')&&isstruct(parameters.sequence)
    sequence=parameters.sequence;
    if isfield(sequence,'name')
        fields(sequence,{'name','count','axis','target'},'parameters.sequence');
        if strcmp(sequence.name,'hahn')&&isfield(sequence,'count')&&sequence.count~=1
            error('cce:parameters','parameters.sequence.count must be one for Hahn echo.');
        end
    else
        fields(sequence,{'type','target','placement','at','axis','angle','operator',...
                         'duration','amplitude','detuning','segments'},'parameters.sequence');
        for event=reshape(sequence,1,[])
            finite=isfield(event,'type')&&strcmp(event.type,'finite');
            if isfield(event,'operator')&&~isempty(event.operator)&&...
               isfield(event,'axis')&&~isempty(event.axis)
                error('cce:parameters','parameters.sequence.axis is unused when operator is supplied.');
            end
            if finite, unused={'angle','operator'}; else, unused={'duration','amplitude','detuning','segments'}; end
            for name=unused
                if isfield(event,name{1})&&~isempty(event.(name{1}))&&...
                   any(event.(name{1})~=0,'all')
                    error('cce:parameters','parameters.sequence.%s is not used by this event type.',name{1});
                end
            end
        end
    end
end

% Check the explicit convergence orchestration vocabulary
if isfield(parameters,'convergence')
    fields(parameters.convergence,{'ladders','atol','rtol','max_work','checks','require_independent'},'parameters.convergence');
    if isfield(parameters.convergence,'ladders')
        fields(parameters.convergence.ladders,{'field','values'},'parameters.convergence.ladders');
    end
    if isfield(parameters.convergence,'checks')
        for k=1:numel(parameters.convergence.checks)
            cce_validate(parameters.convergence.checks{k});
        end
    end
end

% Reject graph scoring auxiliaries outside their implemented scoring modes
if isfield(parameters,'graph')&&isfield(parameters.graph,'score')
    graph=parameters.graph;
    for name={'mismatch','mismatch_scale'}
        if isfield(graph,name{1})&&~strcmp(graph.score,'resonance')
            error('cce:capability','parameters.graph.%s requires score=resonance.',name{1});
        end
    end
    if isfield(graph,'scores')&&~strcmp(graph.score,'supplied')
        error('cce:capability','parameters.graph.scores requires score=supplied.');
    end
end

% Check structured operators and jump records without changing their order
if isfield(parameters,'jumps')
    fields(parameters.jumps,{'support','rate','operator'},'parameters.jumps');
end
if isfield(parameters,'noise')&&isfield(parameters.noise,'terms')
    fields(parameters.noise.terms,{'support','operator'},'parameters.noise.terms');
end

% Keep the numerical backend explicit instead of accepting unused preferences
if isfield(parameters,'numerics')&&isfield(parameters.numerics,'precision')&&...
   ~strcmp(parameters.numerics.precision,'double')
    error('cce:capability','parameters.numerics.precision must be double.');
end

% Bind an explicit propagation preference to the implemented physical method
if isfield(parameters,'numerics')&&isfield(parameters.numerics,'propagator')&&...
   isfield(parameters,'method')
    expected='eig';
    if ismember(parameters.method,{'mecce','megcce'}), expected='stock'; end
    if ~strcmp(parameters.numerics.propagator,expected)
        error('cce:capability','parameters.numerics.propagator must be %s for this method.',expected);
    end
end

% Reject fields whose physical meaning is absent from the requested method
if isfield(parameters,'method')
    conditional=ismember(parameters.method,{'cce','mecce'});
    if ~conditional&&isfield(parameters,'second_order')
        error('cce:capability','parameters.second_order requires cce or mecce.');
    end
    noise=strcmp(parameters.method,'noise');
    if isfield(parameters,'jumps')&&~isempty(parameters.jumps)&&...
       ~ismember(parameters.method,{'mecce','megcce'})
        error('cce:capability','parameters.jumps requires mecce or megcce; Lindblad noise is unsupported.');
    end
    if isfield(parameters,'noise')&&~noise
        error('cce:capability','parameters.noise requires method=noise.');
    end
    if isfield(parameters,'spectrum')&&~noise
        error('cce:capability','parameters.spectrum requires method=noise.');
    end
    if isfield(parameters,'normalized')&&parameters.normalized&&(conditional||noise)
        error('cce:capability','parameters.normalized is only supported for generalized scalar observations.');
    end
end

% Reject unused fit details when fitting is explicitly disabled
if isfield(parameters,'analysis')&&isfield(parameters.analysis,'fit')&&...
   ~parameters.analysis.fit&&numel(fieldnames(parameters.analysis))>1
    error('cce:capability','parameters.analysis fitting options require fit=true.');
end

% Validate explicit raw-retention requests without altering the historical default
if isfield(parameters,'retain_raw')&&...
   (~islogical(parameters.retain_raw)||~isscalar(parameters.retain_raw))
    error('cce:parameters','parameters.retain_raw must be a logical scalar.');
end

% Reject inactive top-level selectors instead of giving them no effect
if isfield(parameters,'operators')&&isfield(parameters,'quantity')&&...
   ~strcmp(parameters.quantity,'observable')
    error('cce:capability','parameters.operators requires quantity=observable.');
end
if isfield(parameters,'meanfield_radius')&&isfield(parameters,'meanfield')&&...
   ~strcmp(parameters.meanfield,'static')
    error('cce:capability','parameters.meanfield_radius requires meanfield=static.');
end

% Reject incompatible bath-state controls for deterministic product states
if isfield(parameters,'bath_state')&&isfield(parameters.bath_state,'mode')
    state=parameters.bath_state; sampled=ismember(state.mode,{'sampled','interlaced'});
    for name={'basis','probabilities','samples','seed'}
        if isfield(state,name{1})&&~sampled
            error('cce:capability','parameters.bath_state.%s requires sampled or interlaced mode.',name{1});
        end
    end
    if isfield(state,'rho')&&~strcmp(state.mode,'explicit_product')
        error('cce:capability','parameters.bath_state.rho requires explicit_product mode.');
    end
end

% Reject unsupported graph selectors in modes that do not consume them
if isfield(parameters,'graph')&&isfield(parameters.graph,'mode')
    graph=parameters.graph;
    for item={'radius','distance';'threshold','coupling';'adjacency','explicit'}.'
        if isfield(graph,item{1})&&~strcmp(graph.mode,item{2})
            error('cce:capability','parameters.graph.%s requires graph.mode=%s.',item{1},item{2});
        end
    end
end

end

% Check a nested field vocabulary before any physical dispatch
function fields(value,allowed,label)
if ~isstruct(value), error('cce:parameters','%s must be a structure.',label); end
unknown=setdiff(fieldnames(value),allowed);
if ~isempty(unknown), error('cce:parameters','Unsupported %s.%s.',label,unknown{1}); end
end

% Require one public parameter structure
function grumble(parameters)
if ~isstruct(parameters)||~isscalar(parameters)
    error('cce:parameters','parameters must be a scalar structure.');
end
end


