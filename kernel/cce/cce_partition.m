% Deterministic physical-spin partitions and bounded exact-patch adaptation.
% Syntax:
%
%                  partition=cce_partition(model,graph,parameters)
%
% Inputs:
%   model, graph - canonical model and bath connectivity.
%   parameters.partition - mode, max_group_size, max_group_dimension, and
%   mode-specific groups or pilot fields. order_semantics is max_groups
%   (default), disjoint, or disjoint_inter; see the source header for pilot
%   score and holdout fields.
%
% Outputs:
%   partition - frozen groups and IDs, group-level graph/order metadata, and
%   pilot/validation diagnostics.
%
% See also cce_clusters, cce_graph.
% Reference: https://doi.org/10.1103/PhysRevB.110.L220302
%
% ilya.kuprov@weizmann.ac.il

function partition=cce_partition(model,graph,parameters)

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate declared bounds and canonicalise physical group membership
grumble(model,graph,parameters);
p=parameters.partition; ids=graph.bath_ids;
[~,positions]=ismember(ids,model.physical_ids); dims=model.multiplicities(positions);
groups=num2cell(ids);
if isfield(parameters,'groups')&&~isempty(parameters.groups)
    groups=cellfun(@sort,parameters.groups,'UniformOutput',false);
end
[~,order]=sort(cellfun(@min,groups)); groups=groups(order);
diagnostics=struct('probes',0,'rounds',0,'reason','',...
                   'pilot_error',[],'heldout_error',[]);
complete=true;

% Merge strongest inter-group edges with physical-ID tie breaking
if strcmp(p.mode,'coupling')
    while true
        edges=zeros(0,4);
        for a=1:numel(groups)
            for b=a+1:numel(groups)
                if dimension([groups{a} groups{b}],ids,dims)>p.max_dimension, continue; end
                [~,ia]=ismember(groups{a},ids); [~,ib]=ismember(groups{b},ids);
                strength=norm(graph.pair_norm(ia,ib),'fro');
                if strength>0, edges(end+1,:)=[-strength min(groups{a}) a b]; end %#ok<AGROW>
            end
        end
        if isempty(edges), break; end
        edges=sortrows(edges,[1 2 4]); a=edges(1,3); b=edges(1,4);
        groups{a}=sort([groups{a} groups{b}]); groups(b)=[];
    end
end

% Assign farthest-first spatial seeds subject to exact Hilbert capacities
if strcmp(p.mode,'spatial')&&~isempty(ids)
    xyz=model.coordinates_A(positions,:); seeds=1;
    for k=2:p.group_count
        distance=inf(numel(ids),1);
        for seed=seeds
            distance=min(distance,sum((xyz-xyz(seed,:)).^2,2));
        end
        distance(seeds)=-inf; [~,next]=max(distance); seeds(end+1)=next; %#ok<AGROW>
    end
    groups=cell(1,p.group_count);
    for k=1:p.group_count, groups{k}=ids(seeds(k)); end
    for k=setdiff(1:numel(ids),seeds,'stable')
        [~,order]=sort(sum((xyz(seeds,:)-xyz(k,:)).^2,2)); placed=false;
        for group=order.'
            if dimension([groups{group} ids(k)],ids,dims)<=p.max_dimension
                groups{group}=sort([groups{group} ids(k)]); placed=true; break;
            end
        end
        if ~placed, error('cce_partition:capacity','Spatial group capacity is insufficient.'); end
    end

    % Accept only strict centroid-objective improvements within the swap budget
    for round=1:p.max_rounds
        improved=false; old=spatial_objective(groups,ids,xyz);
        for a=1:numel(groups)
            for b=a+1:numel(groups)
                for x=groups{a}
                    for y=groups{b}
                        candidate=groups;
                        candidate{a}=sort([setdiff(groups{a},x) y]);
                        candidate{b}=sort([setdiff(groups{b},y) x]);
                        if max(cellfun(@(g)dimension(g,ids,dims),candidate))<=p.max_dimension&&...
                           spatial_objective(candidate,ids,xyz)<old
                            groups=candidate; improved=true; break;
                        end
                    end
                    if improved, break; end
                end
                if improved, break; end
            end
            if improved, break; end
        end
        diagnostics.rounds=round;
        if ~improved, break; end
    end
    if improved, complete=false; diagnostics.reason='spatial_iteration_ceiling'; end
end

% Run optional fixed-pilot adaptation before freezing the production partition
adaptive=strcmp(p.mode,'adaptive')||(isfield(p,'adaptive')&&p.adaptive);
if adaptive&&complete
    [objective,errors,used,valid]=patches(groups,p.pilot_contexts,0);
    diagnostics.probes=used; diagnostics.pilot_error=errors;
    complete=valid&&objective<=1;
    for round=1:p.max_rounds
        if complete||~valid, break; end
        best=objective; accepted=false; selected=groups; selected_errors=errors;
        for a=1:numel(groups)
            for b=a+1:numel(groups)
                if dimension([groups{a} groups{b}],ids,dims)>p.max_dimension, continue; end
                candidate=groups; candidate{a}=sort([groups{a} groups{b}]); candidate(b)=[];
                [value,candidate_errors,used,valid]=patches(candidate,p.pilot_contexts,diagnostics.probes);
                diagnostics.probes=used;
                if ~valid, break; end
                if value<best
                    best=value; selected=candidate; accepted=true; selected_errors=candidate_errors;
                end
            end
            if ~valid, break; end
        end
        diagnostics.rounds=round;
        if ~valid||~accepted, break; end
        groups=selected; objective=best; errors=selected_errors;
        diagnostics.pilot_error=errors; complete=objective<=1;
    end

    % Require all held-out contexts to pass using the frozen pilot-selected groups
    if complete
        [held,errors,used,valid]=patches(groups,p.heldout_contexts,diagnostics.probes);
        diagnostics.heldout_error=errors; diagnostics.probes=used;
        complete=valid&&held<=1;
        if ~complete, diagnostics.reason='heldout_failure_or_probe_ceiling'; end
    else
        diagnostics.reason='patch_tolerance_or_round_probe_dimension_ceiling';
    end
end

% Return canonical groups and expose exhausted validation budgets explicitly
[~,order]=sort(cellfun(@min,groups)); groups=groups(order);
partition=struct('groups',{groups},'dimensions',cellfun(@(g)dimension(g,ids,dims),groups),...
                 'complete',complete,'diagnostics',diagnostics,'mode',p.mode);

function [objective,errors,used,valid]=patches(candidate,contexts,used)

    % Choose the same exact evaluator and group order as production
    objective=0; errors=[]; valid=true;
    q=parameters; q.groups=candidate; q.expansion='standard';
    if strcmp(parameters.expansion,'disjoint'), q.order=1; end
    if strcmp(parameters.expansion,'disjoint_inter'), q.order=2; end
    evaluator=@cce_cluster;
    if ismember(q.method,{'mecce','megcce'}), evaluator=@cce_master; end
    if strcmp(q.method,'noise'), evaluator=@cce_noise; end
    schedule=cce_sequence(model,q);

    % Compare every declared patch with the downward-closed group expansion
    for patch=1:numel(p.patches)
        sites=p.patches{patch};
        if dimension(sites,ids,dims)>p.max_dimension
            valid=false; objective=inf; return;
        end
        local=q; local.bath_spins=sites;
        local.groups=cellfun(@(g)intersect(g,sites),candidate,'UniformOutput',false);
        local.groups=local.groups(~cellfun(@isempty,local.groups));
        local_model=model; local_model.bath_ids=sites;
        local_graph=cce_graph(local_model,local); family=cce_clusters(local_model,local_graph,local);
        empty=struct('id',0,'sites',[]); whole=struct('id',1,'sites',sort(sites));
        for context=1:numel(contexts)
            if used>=p.max_probes, valid=false; objective=inf; return; end
            used=used+1; state=contexts(context);
            baseline=evaluator(model,empty,state,schedule,q);
            exact=evaluator(model,whole,state,schedule,q);
            records=repmat(baseline,1,numel(family));
            for cluster_index=1:numel(family)
                records(cluster_index)=evaluator(model,family(cluster_index),state,schedule,q);
            end

            % Score all requested times without hiding masks or exact-zero tolerances
            reduced=cce_reduce([baseline records],family,local);
            difference=abs(reduced.value-exact.value);
            tolerance=p.patch_atol+p.patch_rtol*abs(exact.value);
            ratio=difference./tolerance; ratio(tolerance==0&difference==0)=0;
            if any(~reduced.valid|~exact.valid|~isfinite(ratio),'all')
                objective=inf;
            else
                objective=max(objective,max(ratio,[],'all'));
            end
            errors(end+1)=max(difference,[],'all'); %#ok<AGROW>
        end
    end

end

end

function value=dimension(sites,ids,dims)

% Measure the exact retained bath Hilbert dimension
[~,index]=ismember(sites,ids); value=prod(dims(index));

end

function objective=spatial_objective(groups,ids,xyz)

% Sum squared Euclidean distances from each physical group's centroid
objective=0;
for k=1:numel(groups)
    [~,index]=ismember(groups{k},ids); coordinates=xyz(index,:);
    objective=objective+sum((coordinates-mean(coordinates,1)).^2,'all');
end

end

function grumble(model,graph,parameters)

if ~isstruct(model)||~isscalar(model)||~all(isfield(model,{'physical_ids','multiplicities'}))||...
   ~isstruct(graph)||~isscalar(graph)||~all(isfield(graph,{'bath_ids','pair_norm'}))||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'partition')||...
   ~isstruct(parameters.partition)||~isscalar(parameters.partition)
    error('Canonical model, graph, and scalar partition parameters are required.');
end
p=parameters.partition; ids=graph.bath_ids;
if ~isnumeric(model.physical_ids)||~isreal(model.physical_ids)||...
   ~isrow(model.physical_ids)||any(~isfinite(model.physical_ids))||...
   ~isnumeric(model.multiplicities)||~isreal(model.multiplicities)||...
   numel(model.multiplicities)~=numel(model.physical_ids)||...
   any(~isfinite(model.multiplicities)|model.multiplicities<1|mod(model.multiplicities,1)~=0)
    error('Physical multiplicities must be finite positive integers.');
end
if ~all(isfield(p,{'mode','max_dimension'}))||~ischar(p.mode)||...
   ~ismember(p.mode,{'user','coupling','spatial','adaptive'})||~integer(p.max_dimension,1)
    error('A recognised partition mode and finite positive integer dimension are required.');
end
if isfield(p,'adaptive')&&(~islogical(p.adaptive)||~isscalar(p.adaptive))
    error('partition.adaptive must be a logical scalar.');
end
if ~isnumeric(ids)||~isreal(ids)||(~isempty(ids)&&~isrow(ids))||...
   any(~isfinite(ids)|ids<1|mod(ids,1)~=0)||numel(unique(ids))~=numel(ids)
    error('Bath IDs must be distinct positive integers in a row.');
end
[found,index]=ismember(ids,model.physical_ids);
if any(~found)||any(model.multiplicities(index)>p.max_dimension)
    error('A bath spin is missing or exceeds the partition capacity.');
end
if ~isnumeric(graph.pair_norm)||~isreal(graph.pair_norm)||...
   ~isequal(size(graph.pair_norm),[numel(ids) numel(ids)])||...
   any(~isfinite(graph.pair_norm)|graph.pair_norm<0,'all')
    error('Physical pair norms must be a finite nonnegative square array.');
end

if isfield(parameters,'groups')&&~isempty(parameters.groups)
    groups=parameters.groups;
    if ~iscell(groups)||~isrow(groups)||any(cellfun(@(g)~isnumeric(g)||...
       ~isreal(g)||~isrow(g)||isempty(g)||any(~isfinite(g)|g<1|mod(g,1)~=0),groups))||...
       ~isequal(sort([groups{:}]),sort(ids))||...
       any(cellfun(@(g)dimension(g,ids,model.multiplicities(index)),groups)>p.max_dimension)
        error('Groups must partition physical bath IDs within capacity.');
    end
end

if strcmp(p.mode,'spatial')
    if ~all(isfield(p,{'group_count','max_rounds'}))||~integer(p.group_count,0)||...
       p.group_count>numel(ids)||(~isempty(ids)&&p.group_count<1)||...
       ~integer(p.max_rounds,1)||~isfield(model,'coordinates_A')||...
       ~isnumeric(model.coordinates_A)||~isreal(model.coordinates_A)||...
       ~isequal(size(model.coordinates_A),[numel(model.physical_ids) 3])||...
       any(~isfinite(model.coordinates_A(index,:)),'all')
        error('Spatial grouping requires finite coordinates, group count, and round ceiling.');
    end
end

if strcmp(p.mode,'adaptive')||(isfield(p,'adaptive')&&p.adaptive)
    needed={'patch_atol','patch_rtol','max_rounds','max_probes','patches',...
            'pilot_contexts','heldout_contexts'};
    if ~all(isfield(p,needed))||~tolerance(p.patch_atol)||~tolerance(p.patch_rtol)||...
       ~integer(p.max_rounds,1)||~integer(p.max_probes,1)||...
       ~iscell(p.patches)||isempty(p.patches)||~isstruct(p.pilot_contexts)||...
       isempty(p.pilot_contexts)||~isstruct(p.heldout_contexts)||isempty(p.heldout_contexts)
        error('Adaptive grouping requires finite bounds, patches, pilot, and held-out contexts.');
    end
    for k=1:numel(p.patches)
        sites=p.patches{k};
        if ~isnumeric(sites)||~isreal(sites)||~isrow(sites)||isempty(sites)||...
           numel(unique(sites))~=numel(sites)||any(~ismember(sites,ids))
            error('Invalid exact patch physical IDs.');
        end
    end
end

end

function valid=integer(value,lower)

% Recognise a finite scalar integer at or above the declared lower bound
valid=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&...
      value>=lower&&mod(value,1)==0;

end

function valid=tolerance(value)

% Recognise a finite nonnegative absolute or relative tolerance
valid=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>=0;

end


