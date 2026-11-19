% Test independent, explicit CCE convergence ladders. Syntax:
%
% [answer,diagnostics]=cce_converge(spin_system,parameters,pulse_sequence,assumptions)
%
% parameters.cce.convergence.ladders has field and values (cell row).
% Supported fields: cce.order, cce.bath_spins, cce.graph.threshold,
% cce.graph.radius, cce.bath_state.samples, times. Each ladder refines the
% baseline independently; no simultaneous ladder changes are allowed.
% parameters.cce.convergence.atol/rtol are required scalar observable tolerances;
% max_work is a positive integer cap on complete clexpn evaluations. A time
% ladder must contain every preceding time exactly (no interpolation).
% parameters.cce.convergence.checks is an optional cell of explicit parameter overrides
% for independent/held-out checks against the baseline. Sampling ladders
% require a held-out seed check. Structural sampling needs explicit bath
% realisations, not an invented generator. Successful output means only
% the requested finite ladder comparisons passed, not a global theorem.
% answer is the unchanged numeric baseline; diagnostics carries the status
% converged, limit_reached, or not_converged; no physics is changed silently.
%
% Parameters:
%
%   spin_system, parameters - physical model and ordinary CCE settings plus
%   cce.convergence.ladders (field and cell-row values), atol, rtol, max_work,
%   and optional checks (cell of explicit parameter overrides).
%   pulse_sequence - registered experiment function handle
%   assumptions - Spinach interaction assumption string
%
% Outputs:
%
%   answer - numeric baseline; diagnostics - convergence status, original
%   input_parameters and convergence_report containing comparisons, work
%   usage and acceptance metadata.
%
% pulse_sequence and assumptions are the explicit clexpn experiment contract.
% Order refinements beyond the finite bath use its full order.
% See also clexpn, cce_diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_converge(spin_system,parameters,pulse_sequence,assumptions)

% Reject shared checkpoints before evaluating any independent refinement
grumble(spin_system,parameters,pulse_sequence,assumptions);
control=parameters.cce.convergence;
baseline=parameters;
baseline.cce=rmfield(baseline.cce,'convergence');
[answer,diagnostics]=clexpn(spin_system,pulse_sequence,baseline,assumptions);
work=1;
model=struct('multiplicities',spin_system.comp.mults,'observed_ids',parameters.observed);
audit=cce_diagnostics(model,diagnostics.raw,diagnostics.reduction,diagnostics.parameters);
capacity=sum(arrayfun(@(x) numel(x.values)-1,control.ladders));
if isfield(control,'checks')
    capacity=capacity+numel(control.checks);
end
report=repmat(struct('field','','step',0,'difference',Inf,'scaled_difference',Inf,...
                     'coverage',false,'numerical_error',Inf,'sampling_error',Inf,...
                     'passed',false),1,capacity);
recorded=0;
all_pass=audit.acceptable;
limited=false;
for k=1:numel(control.ladders)

    % Evaluate successive levels of this independent approximation ladder
    ladder=control.ladders(k);
    passed=true;
    for j=2:numel(ladder.values)
        if work>=control.max_work
            limited=true;
            break;
        end

        % Change exactly one declared approximation using common samples
        trial=baseline;
        pieces=strsplit(ladder.field,'.');
        trial=setfield(trial,pieces{:},ladder.values{j});
        if strcmp(ladder.field,'cce.order')
            trial.cce.order=min(trial.cce.order,numel(diagnostics.bath_ids));
        end
        [current,details]=clexpn(spin_system,pulse_sequence,trial,assumptions);
        work=work+1;
        current_audit=cce_diagnostics(model,details.raw,details.reduction,details.parameters);
        entry=compare(answer,diagnostics,current,details,audit,current_audit,control,ladder.field,j);
        recorded=recorded+1;
        report(recorded)=entry;
        passed=passed&&entry.passed;
    end
    all_pass=all_pass&&passed;
    if limited
        break;
    end

end

% Assess held-out choices separately instead of reusing pilot samples
if isfield(control,'checks')&&~limited
    for k=1:numel(control.checks)
        if work>=control.max_work
            limited=true;
            break;
        end
        trial=baseline;
        override=control.checks{k};
        fields=fieldnames(override);
        for j=1:numel(fields)
            if strcmp(fields{j},'cce')
                controls=fieldnames(override.cce);
                for n=1:numel(controls)
                    trial.cce.(controls{n})=override.cce.(controls{n});
                end
            else
                trial.(fields{j})=override.(fields{j});
            end
        end
        [current,details]=clexpn(spin_system,pulse_sequence,trial,assumptions);
        work=work+1;
        current_audit=cce_diagnostics(model,details.raw,details.reduction,details.parameters);
        entry=compare(answer,diagnostics,current,details,audit,current_audit,control,'independent_check',k);
        recorded=recorded+1;
        report(recorded)=entry;
        all_pass=all_pass&&entry.passed;
    end
end

% Report resource exhaustion distinctly from a failed or accepted ladder
diagnostics.convergence='not_converged';
if limited, diagnostics.convergence='limit_reached';
elseif all_pass
    diagnostics.convergence='converged';
end
diagnostics.input_parameters=parameters;
diagnostics.convergence_report=struct('comparisons',report(1:recorded),'work',work,...
                                      'max_work',control.max_work,'complete_coverage',audit.complete_coverage,...
                                      'atol',control.atol,'rtol',control.rtol,...
                                      'scope','independent finite ladders around the unchanged baseline');

end

% Compare shared times without hiding missing coverage or numerical error
function entry=compare(previous,prior,current,details,audit,current_audit,control,field,step)

% Match the original time mesh and add both numerical error estimates
[shared,ia,ib]=intersect(prior.parameters.times,details.parameters.times,'stable');
coverage=numel(shared)==numel(prior.parameters.times)&&...
         audit.acceptable&&current_audit.acceptable;
if isempty(shared)||size(previous,2)~=size(current,2)
    difference=Inf;
    scaled=Inf;
    numerical=Inf;
    sampling=Inf;
else
    delta=abs(previous(ia,:)-current(ib,:));
    budget=control.atol+control.rtol*max(abs(previous(ia,:)),abs(current(ib,:)));
    numerical=max(prior.reduction.abs_error_bound(ia,:)+details.reduction.abs_error_bound(ib,:),[],'all');
    sampling=max(prior.reduction.sampling_standard_error(ia,:)+...
                 details.reduction.sampling_standard_error(ib,:),...
                 [],'all','includemissing');
    difference=max(delta,[],'all');
    total=delta+prior.reduction.abs_error_bound(ia,:)+details.reduction.abs_error_bound(ib,:);
    ratio=total./budget;
    ratio(total==0&budget==0)=0;
    scaled=max(ratio,[],'all');
    if any(~isfinite(delta),'all')||any(~isfinite(ratio),'all')
        scaled=Inf;
    end
end
entry=struct('field',field,'step',step,'difference',difference,...
             'scaled_difference',scaled,'coverage',coverage,...
             'numerical_error',numerical,'sampling_error',sampling,...
             'passed',coverage&&scaled<=1);

end

% Require explicit independent ladders, tolerances, and a finite work cap
function grumble(spin_system,parameters,pulse_sequence,assumptions)
if ~isa(pulse_sequence,'function_handle')
    error('cce_converge:sequence','Supply an experiment function handle.');
end
if ~(ischar(assumptions)&&isrow(assumptions)&&~isempty(assumptions))&&...
   ~(iscell(assumptions)&&isrow(assumptions)&&~isempty(assumptions)&&ischar(assumptions{1}))
    error('cce_converge:assumptions','Supply explicit stock assume arguments.');
end
if isstruct(parameters)&&isfield(parameters,'checkpoint')
    error('cce_converge:checkpoint','Shared checkpoints are not convergence inputs.');
end
if ~isstruct(spin_system)||~isstruct(parameters)||~isfield(parameters,'cce')||...
   ~isfield(parameters.cce,'convergence')||~isstruct(parameters.cce.convergence)||...
   ~isscalar(parameters.cce.convergence)
    error('cce_converge:input','An explicit convergence structure is required.');
end
control=parameters.cce.convergence;
if ~all(isfield(control,{'ladders','atol','rtol','max_work'}))||...
   ~all(ismember(fieldnames(control),{'ladders','atol','rtol','max_work','checks'}))||...
   ~isstruct(control.ladders)||isempty(control.ladders)||...
   ~all(isfield(control.ladders,{'field','values'}))
    error('cce_converge:input','Explicit ladders, tolerances and max_work are required.');
end
for field={'atol','rtol','max_work'}
    value=control.(field{1});
    if ~isnumeric(value)||~isscalar(value)||~isreal(value)||~isfinite(value)||value<0
        error('cce_converge:input','Convergence budgets must be finite nonnegative scalars.');
    end
end
if control.atol+control.rtol==0||control.max_work<1||control.max_work~=fix(control.max_work)
    error('cce_converge:input','Nonzero tolerance and positive integer max_work required.');
end
if isfield(parameters.cce,'checkpoint')&&~isempty(parameters.cce.checkpoint)
    error('cce_converge:checkpoint','Run independent ladders without a shared checkpoint file.');
end
if numel(unique({control.ladders.field}))~=numel(control.ladders)
    error('cce_converge:input','Each approximation may have only one independent ladder.');
end
allowed={'cce.order','cce.bath_spins','cce.graph.threshold','cce.graph.radius','cce.bath_state.samples','times'};
for k=1:numel(control.ladders)
    ladder=control.ladders(k);
    if ~ischar(ladder.field)||~ismember(ladder.field,allowed)||...
       ~iscell(ladder.values)||numel(ladder.values)<2
        error('cce_converge:input','Each supported ladder needs at least two explicit values.');
    end
    if strcmp(ladder.field,'cce.order')&&...
       ~all(cellfun(@(value)isnumeric(value)&&isreal(value)&&isscalar(value)&&...
                    isfinite(value)&&value>=0&&value==fix(value),ladder.values))
        error('cce_converge:input','Order ladders require finite nonnegative integers.');
    end
    pieces=strsplit(ladder.field,'.');
    initial=parameters;
    for n=1:numel(pieces)
        if ~isstruct(initial)||~isfield(initial,pieces{n})
            error('cce_converge:baseline','A ladder field is absent from the baseline.');
        end
        initial=initial.(pieces{n});
    end
    if ~isequaln(initial,ladder.values{1})
        error('cce_converge:baseline','Every ladder must start at the supplied baseline.');
    end
    if strcmp(ladder.field,'times')
        for j=2:numel(ladder.values)
            if ~all(ismember(ladder.values{j-1},ladder.values{j}))
                error('cce_converge:coverage','Refined times must contain every preceding time.');
            end
        end
    end
end
if isfield(control,'checks')&&(~iscell(control.checks)||...
   ~all(cellfun(@(x) isstruct(x)&&isscalar(x),control.checks)))
    error('cce_converge:input','Independent checks must be a cell of parameter override structures.');
end
if any(strcmp({control.ladders.field},'cce.bath_state.samples'))
    held_out=false;
    if isfield(control,'checks')
        for k=1:numel(control.checks)
            check=control.checks{k};
            held_out=held_out||(isfield(check,'cce')&&isfield(check.cce,'bath_state')&&isfield(check.cce.bath_state,'seed')&&...
                check.cce.bath_state.seed~=parameters.cce.bath_state.seed);
        end
    end
    if ~held_out
        error('cce_converge:sampling','Sampling ladders require an independent seed check.');
    end
end
end


