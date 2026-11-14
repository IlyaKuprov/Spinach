% Deterministic physical-spin partitions and bounded exact-patch adaptation.
% Syntax:
%
%                  partition=cce_partition(model,graph,parameters)
%
% partition.mode: user, coupling, spatial, adaptive. max_dimension is the
% maximum bath-group Hilbert dimension. User groups are parameters.groups.
% Spatial mode also requires group_count and max_rounds. Adaptive mode
% requires patch_atol, patch_rtol, max_rounds, max_probes, patches (cell of
% physical IDs), pilot_contexts and heldout_contexts (cce_states contexts).
% Patches are exact physical spin sets; their total bath dimension must fit
% max_dimension. Production groups are frozen before held-out validation.
% complete=false explicitly reports exhausted ceilings or failed held-out
% validation. Patch success is not a global truncation-error guarantee.
%
% ilya.kuprov@weizmann.ac.il

function partition=cce_partition(model,graph,parameters)

% Validate the physical-ID and bounded grouping contract
grumble(model,graph,parameters);
p=parameters.partition; ids=graph.bath_ids;
[~,positions]=ismember(ids,model.physical_ids); dims=model.multiplicities(positions);
groups=num2cell(ids); diagnostics=struct('probes',0,'rounds',0,'reason','',...
                                      'pilot_error',[],'heldout_error',[]);
complete=true;

% Respect authoritative groups or construct a reproducible bounded partition
switch p.mode
    case 'user'
        groups=parameters.groups;
    case 'coupling'
        if isfield(parameters,'groups')&&~isempty(parameters.groups), groups=parameters.groups; end
        [i,j]=find(triu(graph.pair_norm,1)>0);
        strength=graph.pair_norm(sub2ind(size(graph.pair_norm),i,j));
        [~,order]=sortrows([-strength ids(i).' ids(j).'],[1 2 3]);
        for edge=order.'
            a=find(cellfun(@(g)ismember(ids(i(edge)),g),groups));
            b=find(cellfun(@(g)ismember(ids(j(edge)),g),groups));
            if a~=b&&dimension([groups{a} groups{b}],ids,dims)<=p.max_dimension
                groups{a}=sort([groups{a} groups{b}]); groups(b)=[];
            end
        end
    case 'spatial'
        xyz=model.coordinates_A(positions,:); seeds=1;
        while numel(seeds)<p.group_count
            distance=inf(numel(ids),1);
            for k=seeds, distance=min(distance,sum((xyz-xyz(k,:)).^2,2)); end
            distance(seeds)=-inf; [~,next]=max(distance); seeds(end+1)=next; %#ok<AGROW>
        end
        groups=cell(1,p.group_count);
        for k=1:p.group_count, groups{k}=ids(seeds(k)); end
        for k=setdiff(1:numel(ids),seeds,'stable')
            [~,order]=sort(sum((xyz(seeds,:)-xyz(k,:)).^2,2)); placed=false;
            for j=order.'
                if dimension([groups{j} ids(k)],ids,dims)<=p.max_dimension
                    groups{j}=sort([groups{j} ids(k)]); placed=true; break;
                end
            end
            if ~placed, error('cce_partition:capacity','Spatial group capacity is insufficient.'); end
        end
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
    case 'adaptive'
        if isfield(parameters,'groups')&&~isempty(parameters.groups), groups=parameters.groups; end
        [objective,errors,used,valid]=patches(groups,p.pilot_contexts,0);
        diagnostics.probes=used; diagnostics.pilot_error=errors;
        complete=valid&&objective<=1;
        for round=1:p.max_rounds
            if complete||~valid, break; end
            best=objective; selected=groups; accepted=false;
            for a=1:numel(groups)
                for b=a+1:numel(groups)
                    if dimension([groups{a} groups{b}],ids,dims)>p.max_dimension, continue; end
                    candidate=groups; candidate{a}=sort([groups{a} groups{b}]); candidate(b)=[];
                    [value,candidate_errors,used,valid]=patches(candidate,p.pilot_contexts,diagnostics.probes);
                    diagnostics.probes=used;
                    if ~valid, break; end
                    if value<best
                        best=value; selected=candidate; accepted=true;
                        diagnostics.pilot_error=candidate_errors;
                    end
                end
                if ~valid, break; end
            end
            diagnostics.rounds=round;
            if ~valid||~accepted, break; end
            groups=selected; objective=best; complete=objective<=1;
        end
        if complete
            [held,errors,used,valid]=patches(groups,p.heldout_contexts,diagnostics.probes);
            diagnostics.heldout_error=errors; diagnostics.probes=used;
            complete=valid&&held<=1;
            if ~complete, diagnostics.reason='heldout_failure_or_probe_ceiling'; end
        else
            diagnostics.reason='patch_tolerance_or_round_probe_dimension_ceiling';
        end
end

% Freeze a canonical order and expose rather than hide incomplete convergence
[~,order]=sort(cellfun(@min,groups)); groups=groups(order);
partition=struct('groups',{groups},'complete',complete,'diagnostics',diagnostics,...
                 'mode',p.mode,'dimensions',cellfun(@(g)dimension(g,ids,dims),groups));

    function [objective,errors,used,valid]=patches(candidate,contexts,used)

    % Compare each grouped expansion against an exact local patch trajectory
    objective=0; errors=[]; valid=true;
    q=parameters; q.groups=candidate; q.expansion='standard';
    if strcmp(parameters.expansion,'disjoint'), q.order=1; end
    if strcmp(parameters.expansion,'disjoint_inter'), q.order=2; end
    schedule=cce_sequence(model,q);
    for patch=1:numel(p.patches)
        sites=p.patches{patch};
        if dimension(sites,ids,dims)>p.max_dimension, valid=false; return; end
        local=q; local.groups=cellfun(@(g)intersect(g,sites),candidate,'UniformOutput',false);
        local.groups=local.groups(~cellfun(@isempty,local.groups));
        local.order=min(q.order,numel(local.groups));
        local_model=model; local_model.bath_ids=sort(sites);
        local_graph=cce_graph(local_model,local); family=cce_clusters(local_model,local_graph,local);
        empty=struct('id',0,'sites',[]); whole=struct('id',1,'sites',sort(sites));
        for context=1:numel(contexts)
            if used>=p.max_probes, valid=false; return; end
            used=used+1; state=contexts(context);
            baseline=cce_cluster(model,empty,state,schedule,q);
            exact=cce_cluster(model,whole,state,schedule,q);
            records=repmat(baseline,1,numel(family));
            for k=1:numel(family)
                records(k)=cce_cluster(model,family(k),state,schedule,q);
            end
            reduced=cce_reduce([baseline records],family,local);
            difference=abs(reduced.value-exact.value);
            tolerance=p.patch_atol+p.patch_rtol*abs(exact.value);
            ratio=difference./tolerance;
            ratio(difference==0&tolerance==0)=0;
            if any(~reduced.valid,'all')||any(~isfinite(ratio),'all')
                objective=inf;
            else
                objective=max(objective,max(ratio,[],'all'));
            end
            errors(end+1)=max(difference,[],'all'); %#ok<AGROW>
        end
    end

    end

end

function d=dimension(group,ids,dims)

% Compute a physical Hilbert dimension, never a spin-count surrogate
[~,index]=ismember(group,ids); d=prod(dims(index));

end

function objective=spatial_objective(groups,ids,xyz)

% Sum squared Euclidean distances to each physical group centroid
objective=0;
for k=1:numel(groups)
    [~,index]=ismember(groups{k},ids); coordinates=xyz(index,:);
    objective=objective+sum((coordinates-mean(coordinates,1)).^2,'all');
end

end

function grumble(model,graph,parameters)

% Require declared capacities and valid disjoint physical-ID groups
if ~isstruct(model)||~isstruct(graph)||~isstruct(parameters)||...
   ~isfield(parameters,'partition')||~isstruct(parameters.partition)
    error('Partition parameters are required.');
end
p=parameters.partition; ids=graph.bath_ids;
if ~all(isfield(p,{'mode','max_dimension'}))||...
   ~ismember(p.mode,{'user','coupling','spatial','adaptive'})||...
   ~isscalar(p.max_dimension)||~isfinite(p.max_dimension)||...
   p.max_dimension<1||p.max_dimension~=fix(p.max_dimension)||isempty(ids)
    error('Invalid partition mode or dimension ceiling.');
end
[found,index]=ismember(ids,model.physical_ids);
if ~all(found)||any(model.multiplicities(index)>p.max_dimension), error('Individual site exceeds group capacity.'); end
if isfield(parameters,'groups')&&~isempty(parameters.groups)
    groups=parameters.groups;
    if ~iscell(groups)||any(cellfun(@isempty,groups))||~isequal(sort([groups{:}]),sort(ids))||...
       any(cellfun(@(g)dimension(g,ids,model.multiplicities(index)),groups)>p.max_dimension)
        error('Groups must partition physical bath IDs within capacity.');
    end
elseif strcmp(p.mode,'user'), error('User partition requires groups.'); end
if strcmp(p.mode,'spatial')
    if ~all(isfield(p,{'group_count','max_rounds'}))||...
       p.group_count<1||p.group_count>numel(ids)||p.group_count~=fix(p.group_count)||...
       p.max_rounds<1||p.max_rounds~=fix(p.max_rounds)||any(~isfinite(model.coordinates_A(index,:)),'all')
        error('Spatial grouping requires coordinates, group count and iteration ceiling.');
    end
end
if strcmp(p.mode,'adaptive')
    fields={'patch_atol','patch_rtol','max_rounds','max_probes','patches','pilot_contexts','heldout_contexts'};
    if ~all(isfield(p,fields))||p.patch_atol<0||p.patch_rtol<0||...
       ~isfinite(p.patch_atol+p.patch_rtol)||p.max_rounds<1||p.max_rounds~=fix(p.max_rounds)||...
       p.max_probes<1||p.max_probes~=fix(p.max_probes)||isempty(p.patches)||...
       isempty(p.pilot_contexts)||isempty(p.heldout_contexts)
        error('Adaptive grouping requires bounded pilot, patches and held-out contexts.');
    end
    for k=1:numel(p.patches)
        sites=p.patches{k};
        if isempty(sites)||numel(unique(sites))~=numel(sites)||any(~ismember(sites,ids))
            error('Invalid exact patch physical IDs.');
        end
    end
end

end
