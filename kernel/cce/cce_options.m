% Validate and resolve the native cluster-expansion input contract.
%
% Syntax:
%
% parameters=cce_options(spin_system,parameters,pulse_sequence,assumptions)
%
% Parameters:
%
%   spin_system    - create output; canonical tensors are already in rad/s
%   parameters     - scalar structure with observed (physical spin-ID row),
%                    and times (nonnegative increasing seconds). A transition
%                    selector is required for projected evolution; explicit
%                    generalised-state experiments need not supply one.
%                    A grid name or explicit quadrature replaces orientation;
%                    sum_up selects a weighted answer or unweighted node outputs.
%                    Experiment controls stay at top level; approximation,
%                    cluster selection, ensembles, and reduction live in cce.
%   pulse_sequence - experiment function handle; its detailed preparation,
%                    propagation, and detection checks belong to the experiment
%   assumptions    - explicit stock assume name or cell row of arguments
%
% Outputs:
%
%   parameters - resolved values, never translated old fields. A missing
%                order is returned as [] for the context to resolve.
%
% Defaults: generalised evolution, clusters joined by retained couplings and
% observed-spin mediation, mixed product bath density, no mean field,
% power-counting reduction, and averaging completed expansions.
% Explicit physical populations, quadrupoles,
% polarisation, jump directions/rates, and second-order projection are never
% invented. For closed-system Hilbert Hahn echoes with an energy-level
% transition and ideal midpoint pi pulse, cce.projection='conditional' explicitly drops register-changing
% matrix elements instead of requiring invariant observed manifolds.
% Canonical tensors and declared pump frequencies are not converted.
% noise.descriptors contains operators, spins, and coefficient; its physical
% validation and retained-cluster matrix construction belong to noise_correlation.
% Periodic segments require fixed-order projected free induction, a single
% cycle, and an explicit static exterior product state. Exterior jumps are
% omitted in this approximation; active jumps with split support are refused.
%
% ilya.kuprov@weizmann.ac.il

function parameters=cce_options(spin_system,parameters,pulse_sequence,assumptions)

% Default to a weighted grid sum or to the reference orientation
if isfield(parameters,'grid')
    if ~isfield(parameters,'sum_up')
        parameters.sum_up=true;
    end
    if ischar(parameters.grid)&&isrow(parameters.grid)&&~isempty(parameters.grid)&&...
       ~isfield(parameters,'orientation')
        parameters.grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',parameters.grid),...
                             'alphas','betas','gammas','weights');
    end
elseif ~isfield(parameters,'orientation')
    parameters.orientation=[0 0 0];
end

% Default to generalised evolution of the complete supplied bath
if ~isfield(parameters,'cce')
    parameters.cce=struct();
end
if ~isfield(parameters.cce,'approximation')
    parameters.cce.approximation='generalized';
end
if ~isfield(parameters.cce,'bath_spins')
    parameters.cce.bath_spins=setdiff(1:spin_system.comp.nspins,parameters.observed);
end
if ~isfield(parameters.cce,'order')
    parameters.cce.order=[];
end
if ~isfield(parameters.cce,'mediated')
    parameters.cce.mediated=true;
end

% Default to no mean fields
if ~isfield(parameters.cce,'meanfield')
    parameters.cce.meanfield='none';
end

% Default to power-counting reduction of completed expansions
if ~isfield(parameters.cce,'averaging')
    parameters.cce.averaging='after_expansion';
end
if ~isfield(parameters.cce,'reducer')
    parameters.cce.reducer='powers';
end
if ~isfield(parameters.cce,'backend')
    parameters.cce.backend='auto';
end
if ~isfield(parameters.cce,'expansion')&&isfield(parameters.cce,'groups')
    parameters.cce.expansion='partition';
end
if ~isfield(parameters.cce,'expansion')
    parameters.cce.expansion='standard';
end
if ~isfield(parameters.cce,'retain_raw')
    parameters.cce.retain_raw=false;
end
if ~isfield(parameters.cce,'second_order')
    parameters.cce.second_order=struct('enabled',false);
end

% Default to masked unresolved elements
if ~isfield(parameters.cce,'zero_policy')
    parameters.cce.zero_policy='mask';
end

% Check consistency
grumble(spin_system,parameters,pulse_sequence,assumptions);

end

% Reject unknown fields at a new-contract ownership boundary
function fields(s,allowed,label)
if ~isstruct(s)||~isscalar(s)
    error('cce_options:structure','%s must be a scalar structure.',label);
end
unknown=setdiff(fieldnames(s),allowed);
if ~isempty(unknown)
    error('cce_options:unknown','Unknown %s field: %s.',label,unknown{1});
end
end

% Require an explicitly supported textual choice
function choice(value,allowed,label)
if ~ischar(value)||~isrow(value)||~ismember(value,allowed)
    error('cce_options:choice','Unsupported %s.',label);
end
end

% Validate a physical spin-ID row without changing its order or shape
function spin_ids(ids,count,label,empty_allowed)
if ~isnumeric(ids)||~isreal(ids)||(~isrow(ids)&&~isempty(ids))||...
   (~empty_allowed&&isempty(ids))||any(~isfinite(ids))||any(ids~=fix(ids))||...
   any(ids<1)||any(ids>count)||numel(unique(ids))~=numel(ids)
    error('cce_options:spins','Invalid physical spin-ID row in %s.',label);
end
end

% Validate a true population density, never silently normalise it
function density(rho,dimension)
if ~isnumeric(rho)||~isequal(size(rho),[dimension dimension])||any(~isfinite(rho),'all')
    error('cce_options:density','A density matrix has the wrong dimension or nonfinite entries.');
end
tolerance=dimension*eps*max(1,norm(rho,'fro'));
if norm(rho-rho','fro')>tolerance||abs(trace(rho)-1)>tolerance||...
   min(eig((rho+rho')/2))<-tolerance
    error('cce_options:density','A physical density must be Hermitian, positive and trace one.');
end
end

% Validate a finite nonnegative scalar, optionally allowing positive infinity
function nonnegative(value,label,infinity)
if ~isnumeric(value)||~isreal(value)||~isscalar(value)||isnan(value)||value<0||...
   (~infinity&&~isfinite(value))
    error('cce_options:value','Invalid nonnegative value for %s.',label);
end
end

% Validate an integer option with an explicit lower bound
function integer(value,label,minimum)
if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||...
   value~=fix(value)||value<minimum
    error('cce_options:integer','Invalid integer %s.',label);
end
end

% Validate ownership, required controls and the complete nested vocabulary
function grumble(system,p,experiment,assumptions)
if isfield(p,'sequence')
    error('cce_options:sequence','Use finite_drive, refocusing_train, or hahn_echo, not parameters.sequence.');
end
if isfield(p,'refocusing')&&~isequal(experiment,@hahn_echo)&&~isequal(experiment,@finite_drive)&&~isequal(experiment,@refocusing_train)&&...
   ~isequal(experiment,@cce_compile)&&~isequal(experiment,@cce_generators)
    error('cce_options:refocusing','The refocusing descriptor belongs to hahn_echo.');
end
if isfield(p,'train')&&~isequal(experiment,@refocusing_train)&&~isequal(experiment,@cce_compile)
    error('cce_options:train','Named trains belong to refocusing_train.');
end
if isfield(p,'drive')&&~isequal(experiment,@finite_drive)&&...
   ~isequal(experiment,@cce_compile)&&~isequal(experiment,@cce_generators)
    error('cce_options:drive','The drive descriptor belongs to finite_drive.');
end
fields(p,{'observed','transition','spins','times','nloops','preparation','rho0',...
          'drive','coil','output','normalized','orientation','grid','sum_up','cce','pump',...
          'bath_state','bath_samples','bath_seed','bath_basis','bath_probs',...
          'bath_jumps','noise','spectrum','refocusing','train'},'parameters');
if ~isa(experiment,'function_handle')
    error('cce_options:experiment','pulse_sequence must be a function handle.');
end
if ~(ischar(assumptions)&&isrow(assumptions)&&~isempty(assumptions))&&...
   ~(iscell(assumptions)&&isrow(assumptions)&&~isempty(assumptions)&&ischar(assumptions{1}))
    error('cce_options:assumptions','Supply explicit stock assume arguments.');
end
if ~all(isfield(p,{'observed','times'}))
    error('cce_options:required','observed and times must be supplied.');
end
spin_ids(p.observed,system.comp.nspins,'observed',false);
if ~isnumeric(p.times)||~isreal(p.times)||~isrow(p.times)||isempty(p.times)||...
   any(~isfinite(p.times))||any(p.times<0)||any(diff(p.times)<=0)
    error('cce_options:times','times must be a finite nonnegative strictly increasing row in seconds.');
end
if isfield(p,'grid')
    if isfield(p,'orientation')
        error('cce_options:orientation','Supply either grid or orientation, not both.');
    end
    fields(p.grid,{'alphas','betas','gammas','weights'},'grid');
    grid_check(p.grid);
    weights=double(p.grid.weights);
    if any(weights<0)||~isfinite(sum(weights))||sum(weights)<=0
        error('cce_options:weights','Grid weights must be nonnegative with finite positive sum.');
    end
end
if isfield(p,'sum_up')&&(~(islogical(p.sum_up)||isnumeric(p.sum_up))||...
   ~isreal(p.sum_up)||~isscalar(p.sum_up)||~ismember(p.sum_up,[0 1]))
    error('cce_options:sum_up','sum_up must be a logical scalar or numeric zero or one.');
end
if isfield(p,'orientation')&&(~isnumeric(p.orientation)||~isreal(p.orientation)||...
   ~isequal(size(p.orientation),[1 3])||any(~isfinite(p.orientation)))
    error('cce_options:orientation','orientation must be a finite 1-by-3 radian row.');
end
if isfield(p,'spins')&&(~iscell(p.spins)||~isrow(p.spins)||~all(cellfun(@ischar,p.spins)))
    error('cce_options:channels','spins contains isotope channel names, not physical IDs.');
end
if isfield(p,'nloops')
    integer(p.nloops,'nloops',1);
end
if isfield(p,'normalized')&&(~islogical(p.normalized)||~isscalar(p.normalized))
    error('cce_options:normalized','normalized must be a logical scalar.');
end
if isfield(p,'transition')
    fields(p.transition,{'spins','ms','levels','alpha','beta'},'transition');
    if ~isfield(p.transition,'spins')||~isequal(p.transition.spins,p.observed)
        error('cce_options:transition','transition.spins must equal the ordered observed physical IDs.');
    end
end
if isfield(p,'pump')
    fields(p.pump,{'sites','species','center','halfwidth','units','fraction','allow_empty','boundary'},'pump');
end
if isfield(p,'noise')
    fields(p.noise,{'terms','descriptors','connected','ordering','manifold'},'noise');
    if isfield(p.noise,'terms')
        for k=1:numel(p.noise.terms)
            fields(p.noise.terms(k),{'support','operator'},'noise.terms');
        end
    end
end
if isfield(p,'noise')&&isfield(p.noise,'descriptors')
    for k=1:numel(p.noise.descriptors)
        fields(p.noise.descriptors(k),{'operators','spins','coefficient'},'noise.descriptors');
    end
end
if isfield(p,'spectrum')
    fields(p.spectrum,{'frequency_grid','convention','window','weights','sidedness',...
                       'time_domain','stationary','modulation'},'spectrum');
end
if ~isfield(p,'cce')
    return;
end
fields(p.cce,{'approximation','projection','bath_spins','order','expansion','groups','partition',...
              'mediated','interaction_cutoff','bath_radius',...
              'meanfield','meanfield_radius','interlaced','averaging','reducer',...
              'zero_policy','propagator','backend','second_order','convergence',...
              'checkpoint','retain_raw'},'cce');
contracts={'second_order',{'enabled','gap_tolerance','overlap_tolerance'};...
           'partition',{'mode','max_dimension','group_count','adaptive','patches','patch_atol','patch_rtol','pilot_configs','heldout_configs','max_rounds','max_probes'};...
           'interlaced',{'condition_sites','reference_state','max_states','raw_only'};...
           'checkpoint',{'file','resume'};...
           'convergence',{'ladders','atol','rtol','max_work','checks','require_independent'}};
for k=1:size(contracts,1)
    name=contracts{k,1};
    if isfield(p.cce,name)
        fields(p.cce.(name),contracts{k,2},['cce.' name]);
    end
end
c=p.cce;
ids=p.observed;
bath=c.bath_spins;
nspins=system.comp.nspins;
if isfield(c,'checkpoint')
    if ~isfield(c.checkpoint,'file')||~ischar(c.checkpoint.file)||...
       ~isrow(c.checkpoint.file)||isempty(c.checkpoint.file)||...
       (isfield(c.checkpoint,'resume')&&...
        (~islogical(c.checkpoint.resume)||~isscalar(c.checkpoint.resume)))
        error('cce_options:checkpoint','Checkpoint requests require a nonempty filename and logical resume flag.');
    end
    if isfield(c,'interlaced')||isfield(c,'convergence')
        error('cce_options:checkpoint','Interlacing and independent convergence runs cannot share a checkpoint ledger.');
    end
end
spin_ids(bath,nspins,'parameters.cce.bath_spins',true);
if ~isempty(intersect(ids,bath))
    error('cce_options:ownership','Observed and bath spin IDs must be disjoint.');
end
choice(c.approximation,{'projected','generalized'},'cce.approximation');
choice(c.expansion,{'standard','partition','disjoint','disjoint_inter','periodic_segments'},'cce.expansion');
choice(c.averaging,{'after_expansion','internal_then_expansion'},'cce.averaging');
choice(c.reducer,{'powers','direct'},'cce.reducer');
choice(c.backend,{'auto','direct','stock','analytical_hahn'},'cce.backend');
choice(c.meanfield,{'none','static'},'cce.meanfield');
choice(c.zero_policy,{'mask','error'},'cce.zero_policy');
if ~isempty(c.order)&&(~isnumeric(c.order)||~isreal(c.order)||~isscalar(c.order)||...
   ~isfinite(c.order)||c.order<0||c.order~=fix(c.order)||c.order>numel(bath))
    error('cce_options:order','Order must be an integer within the supplied bath or empty for convergence selection.');
end
for name={'interaction_cutoff','bath_radius','meanfield_radius'}
    if isfield(c,name{1})
        nonnegative(c.(name{1}),['cce.' name{1}],true);
    end
end
if ~islogical(c.retain_raw)||~isscalar(c.retain_raw)
    error('cce_options:retain_raw','retain_raw must be a logical scalar.');
end
if ~isfield(c.second_order,'enabled')||~islogical(c.second_order.enabled)||...
   ~isscalar(c.second_order.enabled)
    error('cce_options:second_order','second_order.enabled must be explicit and logical.');
end
if c.second_order.enabled&&(~strcmp(c.approximation,'projected')||...
   ~all(isfield(c.second_order,{'gap_tolerance','overlap_tolerance'})))
    error('cce_options:second_order','Second-order projection requires projected evolution and explicit gap/overlap criteria.');
end
if isfield(c,'projection')
    choice(c.projection,{'invariant','conditional'},'cce.projection');
    if strcmp(c.projection,'conditional')&&...
       (~strcmp(c.approximation,'projected')||c.second_order.enabled||...
        ~isfield(p,'transition')||~isfield(p.transition,'levels')||isfield(p,'refocusing')||...
        (isfield(p,'bath_jumps')&&~isempty(p.bath_jumps))||...
        (isfield(system,'bas')&&~strcmp(system.bas.formalism,'zeeman-hilb'))||...
        ~(isequal(experiment,@hahn_echo)||isequal(experiment,@cce_compile)))
        error('cce_options:projection','Conditional projection requires a closed-system Hilbert Hahn echo with energy-level selection, an ideal midpoint pulse, and no second-order correction.');
    end
end
for name={'gap_tolerance','overlap_tolerance'}
    if isfield(c.second_order,name{1})
        nonnegative(c.second_order.(name{1}),['cce.second_order.' name{1}],false);
    end
end
if isfield(c,'groups')
    if ~iscell(c.groups)||~isrow(c.groups)||isempty(c.groups)
        error('cce_options:groups','Groups must be a nonempty cell row of physical bath-ID rows.');
    end
    for k=1:numel(c.groups)
        spin_ids(c.groups{k},nspins,'cce.groups',false);
    end
    members=[c.groups{:}];
    if ~isequal(sort(members),sort(bath))
        error('cce_options:groups','Groups must partition the selected bath exactly once.');
    end
end
if strcmp(c.expansion,'periodic_segments')
    if (~isequal(experiment,@free_evolution)&&~isequal(experiment,@cce_compile))||...
       ~strcmp(c.approximation,'projected')||~strcmp(c.meanfield,'static')||...
       ~isfield(p,'bath_state')||~strcmp(c.reducer,'powers')||...
       isempty(c.order)||c.order<1||c.order>=numel(bath)||...
       any(isfield(c,{'groups','partition','convergence','meanfield_radius'}))||...
       (isfield(p,'output')&&~strcmp(p.output,'signal'))
        error('cce_options:periodicSegments',...
              'Periodic segments require a truncated, projected free-induction signal with fixed exterior states.');
    end
elseif ~strcmp(c.expansion,'standard')&&~isfield(c,'groups')&&~isfield(c,'partition')
    error('cce_options:partition','Partitioned expansion needs groups or an explicit partition policy.');
end
if ~islogical(c.mediated)||~isscalar(c.mediated)
    error('cce_options:mediated','cce.mediated must be a logical scalar.');
end
if isfield(p,'bath_state')
    if isfield(p,'bath_samples')
        error('cce_options:state','Supply bath_state density matrices or bath_samples, not both.');
    end
    if ~iscell(p.bath_state)||~isequal(size(p.bath_state),[1 numel(bath)])
        error('cce_options:state','bath_state requires one density matrix per ordered bath ID.');
    end
    for k=1:numel(bath)
        density(p.bath_state{k},system.comp.mults(bath(k)));
    end
end
if isfield(p,'bath_samples')
    if ~isfield(p,'bath_seed')
        error('cce_options:state','Sampled bath states require bath_samples and bath_seed.');
    end
    integer(p.bath_seed,'bath_seed',0);
    if p.bath_seed>2^32-1
        error('cce_options:seed','The seed exceeds the mt19937ar range.');
    end
    integer(p.bath_samples,'bath_samples',1);
    for name={'bath_basis','bath_probs'}
        if isfield(p,name{1})&&(~iscell(p.(name{1}))||~isequal(size(p.(name{1})),[1 numel(bath)]))
            error('cce_options:state','Sampled-state tables must follow the bath IDs.');
        end
    end
    for k=1:numel(bath)
        d=system.comp.mults(bath(k));
        if isfield(p,'bath_basis')
            U=p.bath_basis{k};
            if ~isnumeric(U)||~isequal(size(U),[d d])||any(~isfinite(U),'all')||...
               norm(U'*U-eye(d),'fro')>d*eps*max(1,norm(U,'fro')^2)
                error('cce_options:state','Each sampled basis must be orthonormal.');
            end
        end
        if isfield(p,'bath_probs')
            probabilities=p.bath_probs{k};
            if ~isnumeric(probabilities)||~isreal(probabilities)||...
               ~isequal(size(probabilities),[1 d])||any(~isfinite(probabilities))||...
               any(probabilities<0)||abs(sum(probabilities)-1)>d*eps
                error('cce_options:population','Explicit probabilities must be normalized nonnegative rows.');
            end
        end
    end
elseif any(isfield(p,{'bath_seed','bath_basis','bath_probs'}))
    error('cce_options:state','Sampling fields require bath_samples.');
end
if isfield(c,'interlaced')&&~isfield(p,'bath_samples')
    error('cce_options:state','Interlaced expansion requires sampled bath states.');
end
if isfield(c,'meanfield_radius')&&~strcmp(c.meanfield,'static')
    error('cce_options:meanfield','meanfield_radius requires static mean fields.');
end
if isfield(p,'rho0')
    density(p.rho0,prod(system.comp.mults(ids)));
end
if isfield(p,'rho0')&&isfield(p,'preparation')
    error('cce_options:preparation','Specify a preparation or an explicit rho0, not both.');
end
if isfield(p,'bath_jumps')
    jumps=p.bath_jumps;
    if ~isstruct(jumps)
        error('cce_options:jumps','bath_jumps must be a structure array.');
    end
    for k=1:numel(jumps)
        j=jumps(k);
        fields(j,{'support','operator','rate'},'bath_jumps');
        if ~all(isfield(j,{'support','operator','rate'}))
            error('cce_options:jumps','Every bath jump needs support, operator and rate.');
        end
        spin_ids(j.support,nspins,'bath_jumps.support',false);
        if ~all(ismember(j.support,[ids bath]))
            error('cce_options:jumps','Jump support must belong to selected observed or bath spins.');
        end
        nonnegative(j.rate,'bath_jumps.rate',false);
        if ischar(j.operator)
            choice(j.operator,{'L+','L-','Lx','Ly','Lz'},'bath_jumps.operator');
            if numel(j.support)~=1
                error('cce_options:jumps','A named local jump needs one supporting spin.');
            end
        elseif ~isnumeric(j.operator)||any(~isfinite(j.operator),'all')||...
               ~isequal(size(j.operator),[prod(system.comp.mults(j.support)) prod(system.comp.mults(j.support))])
            error('cce_options:jumps','A jump matrix must match its ordered support.');
        end
    end
end
if isfield(p,'pump')
    if ~all(isfield(p.pump,{'center','halfwidth','units','fraction','allow_empty'}))
        error('cce_options:pump','Pump selection requires center, halfwidth, units, fraction and allow_empty.');
    end
    if ~islogical(p.pump.allow_empty)||~isscalar(p.pump.allow_empty)
        error('cce_options:pump','allow_empty must be a logical scalar.');
    end
    if isfield(p.pump,'sites')
        spin_ids(p.pump.sites,nspins,'pump.sites',p.pump.allow_empty);
        if ~all(ismember(p.pump.sites,bath))
            error('cce_options:pump','Pump sites must belong to the selected bath.');
        end
    end
    if isfield(p.pump,'species')&&(~iscell(p.pump.species)||~isrow(p.pump.species)||...
       ~all(cellfun(@ischar,p.pump.species)))
        error('cce_options:pump','Pump species must be an explicit cell row of isotope names.');
    end
    choice(p.pump.units,{'Hz','rad/s'},'pump.units');
    if isfield(p.pump,'boundary')
        choice(p.pump.boundary,{'open','closed'},'pump.boundary');
    end
    if ~isnumeric(p.pump.center)||~isreal(p.pump.center)||~isscalar(p.pump.center)||~isfinite(p.pump.center)
        error('cce_options:pump','Pump center must be finite in the explicitly declared frequency units.');
    end
    nonnegative(p.pump.halfwidth,'pump.halfwidth',false);
    nonnegative(p.pump.fraction,'pump.fraction',false);
    if p.pump.fraction>1||~islogical(p.pump.allow_empty)||~isscalar(p.pump.allow_empty)
        error('cce_options:pump','Invalid pump fraction or allow_empty flag.');
    end
end
if ~isfield(p,'transition')&&strcmp(c.approximation,'projected')
    error('cce_options:transition','Projected evolution requires an explicit transition selector.');
end
end


