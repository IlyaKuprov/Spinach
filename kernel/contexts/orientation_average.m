% Orientational averaging of a deferred clexpn calculation, without a
% full-bath basis or averaged Hamiltonian. Syntax:
%
%     [answer,diagnostics]=orientation_average(spin_system,task,parameters)
%
% Parameters:
%
%   spin_system     - Spinach system returned by create
%   task            - descriptor returned by context_task(@clexpn,...)
%   parameters.grid - powder grid name, or scalar structure with equally
%                     long real finite alphas, betas, gammas, and weights vectors;
%                     angles are radians, weights nonnegative with positive sum
%
% Outputs:
%
%   answer      - natural experiment output, averaged in the selected order
%   diagnostics - grid, per-orientation diagnostics in orientations, numeric
%                 values in a cell per orientation, and averaged reduction
%
% task.parameters.cce.averaging selects after_expansion (default) or
% internal_then_expansion. The latter averages raw cluster responses and
% their observed-only references before inclusion-exclusion, across both
% orientations and spin configurations; it requires identical cluster families
% at every orientation. Each orientation still constructs its own generators,
% transition selectors, and graph. Zero-weight orientations are omitted.
% Deferred grid weights are normalised. Nested parameters.orientation
% is replaced; all other experiment parameters belong to the deferred task.
% Orientations run serially; cluster/configuration parallelism may own the
% existing pool. Inside an outer worker, nested CCE scheduling is serial.
% A checkpoint file is supported only for one positive-weight orientation.
% Reused spin draws are single statistical samples, not independent draws
% at each orientation; unavailable sampling uncertainties remain NaN.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=orientation_average(spin_system,task,parameters)

% Validate the request before resolving its named quadrature
grumble(spin_system,task,parameters); started=tic;
if ischar(parameters.grid)
    grid=load(fullfile(spin_system.sys.root_dir,'kernel','grids',parameters.grid),...
              'alphas','betas','gammas','weights');
else
    grid=parameters.grid;
end

% Validate the descriptor and resolved quadrature before any nested calculation
grumble(spin_system,task,struct('grid',grid));
fields={'alphas','betas','gammas','weights'};
for n=1:numel(fields), grid.(fields{n})=double(grid.(fields{n})(:)); end
keep=grid.weights>0;
for n=1:numel(fields), grid.(fields{n})=grid.(fields{n})(keep); end
grid.weights=grid.weights/sum(grid.weights); count=numel(grid.weights);
if count>1&&isfield(task.parameters,'cce')&&isfield(task.parameters.cce,'checkpoint')
    error('orientation_average:checkpoint',...
          'Checkpoint files require exactly one positive-weight orientation.');
end

% Keep parallel ownership with an outer worker rather than starting a nested pool
if ~isempty(getCurrentTask())
    [task.parameters,~]=cce_options(spin_system,task.parameters,...
                                   task.pulse_sequence,task.assumptions);
    task.parameters.cce.parallel.mode='serial';
end

% Run independently rotated physical systems without a full-bath basis
orientations=cell(1,count); values=cell(1,count);
for n=1:count
    nested=task.parameters;
    nested.orientation=[grid.alphas(n) grid.betas(n) grid.gammas(n)];
    [values{n},orientations{n}]=task.context(spin_system,...
                                           task.pulse_sequence,nested,task.assumptions);
end

% Select the prescribed ordering of orientational and cluster reductions
parameters=orientations{1}.parameters;
reduction=orientations{1}.reduction;
if strcmp(parameters.cce.averaging,'internal_then_expansion')
    raw={};
    for n=1:count
        if ~isequal(orientations{n}.clusters,orientations{1}.clusters)
            error('orientation_average:clusters',...
                  'Internal averaging requires identical cluster families at every orientation.');
        end
        records=orientations{n}.raw;
        for s=1:numel(records)
            for k=1:numel(records{s})
                records{s}(k).weight=grid.weights(n)*records{s}(k).weight;
            end
        end
        raw=[raw records]; %#ok<AGROW>
    end
    reduction=cce_reduce(raw,orientations{1}.clusters,parameters);

    % Quadrature nodes sharing a spin draw are not independent statistical samples
    reduction.sampling_standard_error=orientations{1}.reduction.sampling_standard_error;
    samples=numel(orientations{1}.raw);
    complete=find(arrayfun(@(item)isequal(sort(item.sites),...
                                        sort(parameters.cce.bath_spins)),orientations{1}.clusters),1);
    if strcmp(parameters.cce.bath_state.mode,'sampled')&&samples>1&&...
       ~isempty(complete)&&~(isfield(parameters,'noise')&&parameters.noise.connected)
        state_values=zeros([size(reduction.value) samples]);
        for n=1:count
            for k=1:samples
                state_values(:,:,k)=state_values(:,:,k)+grid.weights(n)*...
                                      orientations{n}.raw{k}(complete+1).value;
            end
        end
        reduction.sampling_standard_error=std(state_values,0,3)/sqrt(samples);
    end
else

    % Accumulate completed expansions and conservative absolute error bounds
    fields={'value','by_order','abs_error_bound','by_order_error'};
    for k=1:numel(fields), reduction.(fields{k})=0*reduction.(fields{k}); end
    for n=1:count
        current=orientations{n}.reduction;
        for k=1:numel(fields)
            name=fields{k};
            reduction.(name)=reduction.(name)+grid.weights(n)*current.(name);
        end
        reduction.valid=reduction.valid&current.valid;
        reduction.by_order_valid=reduction.by_order_valid&current.by_order_valid;
        failed=~current.by_order_valid;
        reduction.by_order_reason(failed)=current.by_order_reason(failed);
    end
    reduction.reason=reduction.by_order_reason(:,:,end);

    % Reuse the same sampled spin states across orientations before estimating error
    samples=reduction.samples;
    for k=1:numel(samples)
        for j=1:numel(fields), samples{k}.(fields{j})=0*samples{k}.(fields{j}); end
        for n=1:count
            current=orientations{n}.reduction;
            if ~isequal(current.weights,reduction.weights)||numel(current.samples)~=numel(samples)
                error('orientation_average:states','Orientation sampling requires a common spin ensemble.');
            end
            for j=1:numel(fields)
                name=fields{j};
                samples{k}.(name)=samples{k}.(name)+grid.weights(n)*current.samples{k}.(name);
            end
            samples{k}.valid=samples{k}.valid&current.samples{k}.valid;
            samples{k}.by_order_valid=samples{k}.by_order_valid&current.samples{k}.by_order_valid;
            failed=~current.samples{k}.by_order_valid;
            samples{k}.by_order_reason(failed)=current.samples{k}.by_order_reason(failed);
        end
    end
    for k=1:numel(samples), samples{k}.reason=samples{k}.by_order_reason(:,:,end); end
    reduction.samples=samples;
    if numel(samples)>1&&all(abs(reduction.weights-1/numel(samples))<=eps)
        sampled=cellfun(@(x)x.value,samples,'UniformOutput',false);
        reduction.sampling_standard_error=std(cat(3,sampled{:}),0,3)/sqrt(numel(samples));
    end
end

% Return the numeric answer separately from orientation-resolved evidence
answer=reduction.value; diagnostics.grid=grid;
diagnostics.orientations=orientations; diagnostics.values=values;
diagnostics.parameters=parameters; diagnostics.reduction=reduction;
diagnostics.elapsed_seconds=toc(started);

end

% Validate all descriptor fields and the grid source
function grumble(spin_system,task,parameters)
if ~isstruct(spin_system)||~isscalar(spin_system)||~isfield(spin_system,'comp')
    error('orientation_average:system','Supply a Spinach system returned by create.');
end
fields={'context','pulse_sequence','parameters','assumptions'};
if ~isstruct(task)||~isscalar(task)||~all(isfield(task,fields))||...
   numel(fieldnames(task))~=numel(fields)
    error('orientation_average:task','Supply a context_task descriptor.');
end
if ~isa(task.context,'function_handle')||~isequal(task.context,@clexpn)||...
   ~isa(task.pulse_sequence,'function_handle')||...
   ~isstruct(task.parameters)||~isscalar(task.parameters)||...
   ~ischar(task.assumptions)||~isrow(task.assumptions)||isempty(task.assumptions)
    error('orientation_average:task','Supply a valid deferred clexpn experiment.');
end
if isfield(task.parameters,'sequence')
    error('context_task:sequence','Use finite_drive, refocusing_train, or hahn_echo, not parameters.sequence.');
end
cce_capabilities(task.pulse_sequence);
if ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'grid')||...
   ~(isstruct(parameters.grid)||(ischar(parameters.grid)&&...
     isrow(parameters.grid)&&~isempty(parameters.grid)))
    error('orientation_average:grid','Supply a grid name or explicit quadrature structure.');
end
if ischar(parameters.grid), return; end
grid=parameters.grid;
fields={'alphas','betas','gammas','weights'};
if ~isstruct(grid)||~isscalar(grid)||~all(isfield(grid,fields))||...
   numel(fieldnames(grid))~=numel(fields)
    error('orientation_average:grid','Grid fields must be alphas, betas, gammas and weights.');
end
count=numel(grid.weights);
for n=1:numel(fields)
    value=grid.(fields{n});
    if ~isnumeric(value)||~isreal(value)||~isvector(value)||...
       numel(value)~=count||any(~isfinite(value))
        error('orientation_average:grid','Grid fields must be equally long finite real vectors.');
    end
end
if any(grid.weights<0)||~isfinite(sum(grid.weights))||sum(grid.weights)<=0
    error('orientation_average:weights','Quadrature weights must be nonnegative with positive sum.');
end
end


