% Test independent, explicit CCE convergence ladders. Syntax:
%
%                  result=cce_converge(spin_system,parameters)
%
% convergence.ladders is a struct array with field and values (cell row).
% Supported fields: order, bath_spins, graph.threshold, graph.radius,
% bath_state.samples, times. Each ladder refines one choice from the
% baseline independently; no simultaneous parameter changes are allowed.
% convergence.atol/rtol are required scalar observable tolerances;
% max_work is a positive integer cap on complete cce evaluations. A time
% ladder must contain every preceding time exactly (no interpolation).
% convergence.checks is an optional cell of explicit parameter overrides
% for independent/held-out checks against the baseline. Sampling ladders
% require a held-out seed check. Structural sampling needs explicit bath
% realisations, not an invented generator. Successful output means only
% the requested finite ladder comparisons passed, not a global theorem.
% result is the baseline result with convergence_report and status
% converged, limit_reached or not_converged; no physics is changed silently.
%
% Inputs:
%   spin_system, parameters - physical model and ordinary CCE settings plus
%   convergence.ladders (field and cell-row values), atol, rtol, max_work,
%   and optional checks (cell of complete parameter structures).
%
% Outputs:
%   result - baseline CCE result with convergence status, original
%   input_parameters and convergence_report containing comparisons, work
%   usage and acceptance metadata.
%
% ilya.kuprov@weizmann.ac.il

function result=cce_converge(spin_system,parameters)

% Validate all ladder fields and budgets before running the baseline
grumble(spin_system,parameters);
control=parameters.convergence; baseline=rmfield(parameters,'convergence');
result=cce(spin_system,baseline); work=1;
capacity=sum(arrayfun(@(x) numel(x.values)-1,control.ladders));
if isfield(control,'checks'), capacity=capacity+numel(control.checks); end
report=repmat(struct('field','','step',0,'difference',Inf,'scaled_difference',Inf,...
                     'coverage',false,'numerical_error',Inf,'sampling_error',Inf,...
                     'passed',false),1,capacity); recorded=0;
all_pass=result.diagnostics.acceptable; limited=false;
for k=1:numel(control.ladders)
    ladder=control.ladders(k); previous=result;
    passed=false;
    for j=2:numel(ladder.values)
        if work>=control.max_work, limited=true; break; end

        % Change exactly one declared approximation using common samples
        trial=baseline; pieces=strsplit(ladder.field,'.');
        if isscalar(pieces), trial.(pieces{1})=ladder.values{j};
        else, trial.(pieces{1}).(pieces{2})=ladder.values{j}; end
        current=cce(spin_system,trial); work=work+1;
        entry=compare(previous,current,control,ladder.field,j);
        recorded=recorded+1; report(recorded)=entry; passed=entry.passed; previous=current;
    end
    all_pass=all_pass&&passed;
    if limited, break; end
end

% Assess held-out choices separately instead of reusing pilot samples
if isfield(control,'checks')&&~limited
    for k=1:numel(control.checks)
        if work>=control.max_work, limited=true; break; end
        trial=baseline; override=control.checks{k}; fields=fieldnames(override);
        for j=1:numel(fields), trial.(fields{j})=override.(fields{j}); end
        current=cce(spin_system,trial); work=work+1;
        entry=compare(result,current,control,'independent_check',k);
        recorded=recorded+1; report(recorded)=entry; all_pass=all_pass&&entry.passed;
    end
end

% Report resource exhaustion distinctly from a failed or accepted ladder
result.convergence='not_converged';
if limited, result.convergence='limit_reached';
elseif all_pass, result.convergence='converged'; end
result.input_parameters=parameters;
result.convergence_report=struct('comparisons',report(1:recorded),'work',work,...
    'max_work',control.max_work,'complete_coverage',result.diagnostics.complete_coverage,...
    'atol',control.atol,'rtol',control.rtol,...
    'scope','independent finite ladders around the unchanged baseline');

end

% Compare shared times without hiding missing coverage or numerical error
function entry=compare(previous,current,control,field,step)
[shared,ia,ib]=intersect(previous.times,current.times,'stable');
coverage=numel(shared)==numel(previous.times)&&...
         previous.diagnostics.acceptable&&current.diagnostics.acceptable;
if isempty(shared)||size(previous.value,2)~=size(current.value,2)
    difference=Inf; scaled=Inf; numerical=Inf; sampling=Inf;
else
    delta=abs(previous.value(ia,:)-current.value(ib,:));
    budget=control.atol+control.rtol*max(abs(previous.value(ia,:)),abs(current.value(ib,:)));
    numerical=max(previous.abs_error_bound(ia,:)+current.abs_error_bound(ib,:),[],'all');
    sampling=max(previous.sampling_standard_error(ia,:)+current.sampling_standard_error(ib,:),[],'all');
    difference=max(delta,[],'all');
    total=delta+previous.abs_error_bound(ia,:)+current.abs_error_bound(ib,:);
    ratio=total./budget; ratio(total==0&budget==0)=0;
    scaled=max(ratio,[],'all');
    if any(~isfinite(delta),'all')||any(~isfinite(ratio),'all'), scaled=Inf; end
end
entry=struct('field',field,'step',step,'difference',difference,...
             'scaled_difference',scaled,'coverage',coverage,...
             'numerical_error',numerical,'sampling_error',sampling,...
             'passed',coverage&&scaled<=1);
end

% Require explicit independent ladders, tolerances and a finite work cap
function grumble(spin_system,parameters)
if ~isstruct(spin_system)||~isstruct(parameters)||~isfield(parameters,'convergence')||...
   ~isstruct(parameters.convergence)||~isscalar(parameters.convergence)
    error('cce_converge:input','An explicit convergence structure is required.');
end
control=parameters.convergence;
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
if isfield(parameters,'checkpoint')&&~isempty(parameters.checkpoint)
    error('cce_converge:checkpoint','Run independent ladders without a shared checkpoint file.');
end
if numel(unique({control.ladders.field}))~=numel(control.ladders)
    error('cce_converge:input','Each approximation may have only one independent ladder.');
end
allowed={'order','bath_spins','graph.threshold','graph.radius','bath_state.samples','times'};
for k=1:numel(control.ladders)
    ladder=control.ladders(k);
    if ~ischar(ladder.field)||~ismember(ladder.field,allowed)||...
       ~iscell(ladder.values)||numel(ladder.values)<2
        error('cce_converge:input','Each supported ladder needs at least two explicit values.');
    end
    pieces=strsplit(ladder.field,'.');
    if ~isfield(parameters,pieces{1})
        error('cce_converge:baseline','A ladder field is absent from the baseline.');
    end
    initial=parameters.(pieces{1});
    if ~isscalar(pieces)
        if ~isstruct(initial)||~isfield(initial,pieces{2})
            error('cce_converge:baseline','A nested ladder field is absent from the baseline.');
        end
        initial=initial.(pieces{2});
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
if any(strcmp({control.ladders.field},'bath_state.samples'))
    held_out=false;
    if isfield(control,'checks')
        for k=1:numel(control.checks)
            check=control.checks{k};
            held_out=held_out||(isfield(check,'bath_state')&&isfield(check.bath_state,'seed')&&...
                check.bath_state.seed~=parameters.bath_state.seed);
        end
    end
    if ~held_out, error('cce_converge:sampling','Sampling ladders require an independent seed check.'); end
end
end

