% Deterministic physical-spin partitions and bounded exact-patch adaptation.
% Syntax:
%
%                  partition=cce_partition(model,graph,parameters)
%
% Parameters:
%
%   model, graph - canonical model and bath connectivity.
%   parameters.cce.partition - native user/coupling/spatial/adaptive policy
%   (maximum group Hilbert dimension). Spatial mode also requires
%   group_count and max_rounds. Adaptive mode or optional adaptive=true requires patches,
%   patch_atol, patch_rtol, pilot_configs, heldout_configs, max_rounds,
%   and max_probes. Explicit groups, expansion, and order are CCE options,
%   not partition fields.
%
% Outputs:
%
%   partition - frozen groups and IDs, group-level graph/order metadata, and
%   pilot/validation diagnostics.
%
% See also cce_clusters, cce_graph.
% Reference: https://doi.org/10.1103/PhysRevB.110.L220302
%
% ilya.kuprov@weizmann.ac.il

function partition=cce_partition(model,graph,parameters)

% Validate native partition controls before constructing internal selections
grumble(model,graph,parameters);

% Keep experiment parameters intact while resolving the partition controls
native_parameters=parameters;
options=parameters.cce;
parameters=struct('partition',options.partition,'graph',options.graph,...
                  'order',options.order,'numerics',options.numerics,...
                  'expansion',options.expansion,'method','gcce');
if isfield(parameters.numerics,'propagator')
    parameters.numerics=rmfield(parameters.numerics,'propagator');
end
if isfield(options,'groups')
    parameters.groups=options.groups;
end
if strcmp(options.approximation,'projected')
    parameters.method='cce';
end
if isfield(native_parameters,'noise')
    parameters.method='noise';
    parameters.noise=native_parameters.noise;
    parameters.bath_state=options.bath_state;
end

% Canonicalise physical group membership
p=parameters.partition;
ids=graph.bath_ids;
[~,positions]=ismember(ids,model.physical_ids);
dims=model.multiplicities(positions);
groups=num2cell(ids);
if isfield(parameters,'groups')&&~isempty(parameters.groups)
    groups=cellfun(@sort,parameters.groups,'UniformOutput',false);
end
[~,order]=sort(cellfun(@min,groups));
groups=groups(order);
diagnostics=struct('probes',0,'rounds',0,'reason','',...
                   'pilot_error',[],'heldout_error',[]);
complete=true;

% Merge strongest inter-group edges with physical-ID tie breaking
if strcmp(p.mode,'coupling')
    while true
        edges=zeros(0,4);
        [~,sites]=ismember([groups{:}],ids);
        ng=numel(groups);
        vertices=cell(size(groups));
        for k=1:ng
            vertices{k}=repmat(k,size(groups{k}));
        end
        group=sparse(sites,[vertices{:}],1,numel(ids),ng);
        [left,right]=find(triu(group.'*spones(graph.pair_norm)*group,1));
        for k=1:numel(left)
            a=left(k);
            b=right(k);
            if dimension([groups{a} groups{b}],ids,dims)<=p.max_dimension
                [~,ia]=ismember(groups{a},ids);
                [~,ib]=ismember(groups{b},ids);
                strength=norm(graph.pair_norm(ia,ib),'fro');
                edges(end+1,:)=[-strength min(groups{a}) a b]; %#ok<AGROW>
            end
        end
        if isempty(edges)
            break;
        end
        edges=sortrows(edges,[1 2 4]);
        a=edges(1,3);
        b=edges(1,4);
        groups{a}=sort([groups{a} groups{b}]);
        groups(b)=[];
    end
end

% Assign farthest-first spatial seeds subject to exact Hilbert capacities
if strcmp(p.mode,'spatial')&&~isempty(ids)
    xyz=model.coordinates_A(positions,:);
    seeds=1;
    for k=2:p.group_count
        distance=inf(numel(ids),1);
        for j=seeds
            distance=min(distance,sum((xyz-xyz(j,:)).^2,2));
        end
        distance(seeds)=-inf;
        [~,next]=max(distance);
        seeds(end+1)=next; %#ok<AGROW>
    end
    groups=cell(1,p.group_count);
    for k=1:p.group_count
        groups{k}=ids(seeds(k));
    end
    for k=setdiff(1:numel(ids),seeds,'stable')
        [~,order]=sort(sum((xyz(seeds,:)-xyz(k,:)).^2,2));
        placed=false;
        for g=order.'
            if dimension([groups{g} ids(k)],ids,dims)<=p.max_dimension
                groups{g}=sort([groups{g} ids(k)]);
                placed=true;
                break;
            end
        end
        if ~placed
            error('cce_partition:capacity','Spatial group capacity is insufficient.');
        end
    end

    % Accept only strict centroid-objective improvements within the swap budget
    for r=1:p.max_rounds
        improved=false;
        old=spatial_objective(groups,ids,xyz);
        for a=1:numel(groups)
            for b=a+1:numel(groups)
                for x=groups{a}
                    for y=groups{b}
                        candidate=groups;
                        candidate{a}=sort([setdiff(groups{a},x) y]);
                        candidate{b}=sort([setdiff(groups{b},y) x]);
                        if max(cellfun(@(g)dimension(g,ids,dims),candidate))<=p.max_dimension&&...
                           spatial_objective(candidate,ids,xyz)<old
                            groups=candidate;
                            improved=true;
                            break;
                        end
                    end
                    if improved
                        break;
                    end
                end
                if improved
                    break;
                end
            end
            if improved
                break;
            end
        end
        diagnostics.rounds=r;
        if ~improved
            break;
        end
    end
    if improved
        complete=false;
        diagnostics.reason='spatial_iteration_ceiling';
    end
end

% Run optional fixed-pilot adaptation before freezing the production partition
adaptive=strcmp(p.mode,'adaptive')||(isfield(p,'adaptive')&&p.adaptive);
if adaptive&&complete

    % Both cohorts remain resident throughout pilot selection and validation
    pilot=p.pilot_configs;
    heldout=p.heldout_configs;
    metadata=whos('pilot','heldout');
    count=min(numel(p.patches)*(numel(pilot)+numel(heldout)),p.max_probes);
    bytes=sum([metadata.bytes])+3*8*count;
    if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
        error('cce_partition:memory',...
              'Pilot and held-out configurations exceed numerics.max_memory_bytes.');
    end
    [objective,errors,used,valid]=patches(groups,p.pilot_configs,0);
    diagnostics.probes=used;
    diagnostics.pilot_error=errors;
    complete=valid&&objective<=1;
    for r=1:p.max_rounds
        if complete||~valid
            break;
        end
        best=objective;
        accepted=false;
        selected=groups;
        selected_errors=errors;
        for a=1:numel(groups)
            for b=a+1:numel(groups)
                if dimension([groups{a} groups{b}],ids,dims)>p.max_dimension
                    continue;
                end
                candidate=groups;
                candidate{a}=sort([groups{a} groups{b}]);
                candidate(b)=[];
                [value,candidate_errors,used,valid]=patches(candidate,p.pilot_configs,diagnostics.probes);
                diagnostics.probes=used;
                if ~valid
                    break;
                end

                % Allow bounded progress through unresolved intermediate reductions
                if value<best||(isinf(value)&&isinf(best)&&~accepted)
                    best=value;
                    selected=candidate;
                    accepted=true;
                    selected_errors=candidate_errors;
                end
            end
            if ~valid
                break;
            end
        end
        diagnostics.rounds=r;
        if ~valid||~accepted
            break;
        end
        groups=selected;
        objective=best;
        errors=selected_errors;
        diagnostics.pilot_error=errors;
        complete=objective<=1;
    end

    % Require all held-out spin_configs to pass using the frozen pilot-selected groups
    if complete
        [held,errors,used,valid]=patches(groups,p.heldout_configs,diagnostics.probes);
        diagnostics.heldout_error=errors;
        diagnostics.probes=used;
        complete=valid&&held<=1;
        if ~complete
            diagnostics.reason='heldout_failure_or_probe_ceiling';
        end
    else
        diagnostics.reason='patch_tolerance_or_round_probe_dimension_ceiling';
    end

end

% Return canonical groups and expose exhausted validation budgets explicitly
[~,order]=sort(cellfun(@min,groups));
groups=groups(order);
partition=struct('groups',{groups},'dimensions',cellfun(@(g)dimension(g,ids,dims),groups),...
                 'complete',complete,'diagnostics',diagnostics,'mode',p.mode);

% Evaluate exact patches against the production group expansion
function [objective,errors,used,valid]=patches(candidate,spin_configs,used)

    % Choose the same exact evaluator and group order as production
    objective=0;
    errors=[];
    valid=true;
    q=parameters;
    q.groups=candidate;
    q.expansion='standard';

    % Reserve checkpoint ownership for the production scheduler
    if isfield(q,'checkpoint')
        q=rmfield(q,'checkpoint');
    end
    if strcmp(parameters.expansion,'disjoint')
        q.order=1;
    end
    if strcmp(parameters.expansion,'disjoint_inter')
        q.order=2;
    end
    contract=cce_capabilities(model.pulse_sequence);
    evaluator=contract.evaluator;
    schedule=struct('pulse_sequence',model.pulse_sequence,'contract',contract);
    evaluation=native_parameters;
    evaluation.cce.groups=q.groups;
    evaluation.cce.order=q.order;
    evaluation.cce.expansion=q.expansion;
    if isfield(evaluation.cce,'checkpoint')
        evaluation.cce=rmfield(evaluation.cce,'checkpoint');
    end

    % Compare every declared patch with the downward-closed group expansion
    for n=1:numel(p.patches)
        sites=p.patches{n};
        if dimension(sites,ids,dims)>p.max_dimension
            valid=false;
            objective=inf;
            return;
        end
        local=q;
        local.bath_spins=sites;
        local.groups=cellfun(@(g)intersect(g,sites),candidate,'UniformOutput',false);
        local.groups=local.groups(~cellfun(@isempty,local.groups));

        % Restrict bath-indexed graph data in the declared physical patch order
        [~,indices]=ismember(sites,model.bath_ids);
        for field={'adjacency','scores','mismatch'}
            if isfield(local.graph,field{1})
                matrix=local.graph.(field{1});
                local.graph.(field{1})=matrix(indices,indices);
            end
        end
        local_model=model;
        local_model.bath_ids=sites;
        local_graph=cce_graph(local_model,local);
        family=cce_clusters(local_model,local_graph,local);
        dc=prod(model.multiplicities(model.observed_ids));
        components=1;
        if isfield(native_parameters,'output')&&strcmp(native_parameters.output,'density')
            components=dc^2;
        elseif isfield(native_parameters,'coil')&&iscell(native_parameters.coil)
            components=numel(native_parameters.coil);
        end
        times=native_parameters.times;
        bytes=3*(numel(family)+2)*...
              (numel(times)*(25*components+16*dc^2)+1024);
        if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
            error('cce_partition:memory',...
                  'Patch cluster records exceed numerics.max_memory_bytes.');
        end
        empty=struct('id',0,'sites',[]);
        whole=struct('id',1,'sites',sort(sites));
        full=find(arrayfun(@(entry)isequal(sort(entry.sites),sort(sites)),family),1);
        for s=1:numel(spin_configs)
            if used>=p.max_probes
                valid=false;
                objective=inf;
                return;
            end
            used=used+1;
            state=spin_configs(s);
            state.weight=1;
            baseline=evaluator(model,empty,state,schedule,evaluation);
            exact=evaluator(model,whole,state,schedule,evaluation);
            records=repmat(baseline,1,numel(family));
            for c=1:numel(family)
                records(c)=evaluator(model,family(c),state,schedule,evaluation);
            end

            % Score all times using exact complete patches without auxiliary denominators
            if ~isempty(full)
                reduced=records(full);
            elseif isempty(sites)
                reduced=baseline;
            else
                request=evaluation;
                request.cce.groups=local.groups;
                request.cce.bath_spins=sites;
                reduced=cce_reduce({[baseline records]},family,request);
            end
            difference=abs(reduced.value-exact.value);
            tolerance=p.patch_atol+p.patch_rtol*abs(exact.value);
            ratio=difference./tolerance;
            ratio(tolerance==0&difference==0)=0;
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

% Measure the exact retained bath Hilbert dimension
function value=dimension(sites,ids,dims)
[~,index]=ismember(sites,ids);
value=prod(dims(index));

end

% Sum squared Euclidean distances from each physical group's centroid
function objective=spatial_objective(groups,ids,xyz)
objective=0;
for k=1:numel(groups)
    [~,index]=ismember(groups{k},ids);
    coordinates=xyz(index,:);
    objective=objective+sum((coordinates-mean(coordinates,1)).^2,'all');
end

end

% Recognise a finite scalar integer at or above the declared lower bound
function valid=integer(value,lower)
valid=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&...
      value>=lower&&mod(value,1)==0;

end

% Recognise a finite nonnegative absolute or relative tolerance
function valid=tolerance(value)
valid=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>=0;

end

% Validate declared group membership and partition bounds
function grumble(model,graph,parameters)
if ~isstruct(model)||~isscalar(model)||~all(isfield(model,{'physical_ids','multiplicities'}))||...
   ~isstruct(graph)||~isscalar(graph)||~all(isfield(graph,{'bath_ids','pair_norm'}))||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'cce')||...
   ~isstruct(parameters.cce)||~isscalar(parameters.cce)||...
   ~all(isfield(parameters.cce,{'partition','graph','order','numerics','expansion','approximation'}))||...
   ~isstruct(parameters.cce.partition)||~isscalar(parameters.cce.partition)
    error('Canonical model, graph, and scalar partition parameters are required.');
end
parameters=parameters.cce;
p=parameters.partition;
ids=graph.bath_ids;
if ~isstruct(parameters.numerics)||~isscalar(parameters.numerics)||...
   ~isfield(parameters.numerics,'max_memory_bytes')||...
   ~isnumeric(parameters.numerics.max_memory_bytes)||...
   ~isreal(parameters.numerics.max_memory_bytes)||...
   ~isscalar(parameters.numerics.max_memory_bytes)||...
   ~isfinite(parameters.numerics.max_memory_bytes)||parameters.numerics.max_memory_bytes<=0
    error('A finite positive partition memory budget is required.');
end
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
   any(~isfinite(nonzeros(graph.pair_norm)))||any(graph.pair_norm<0,'all')
    error('Physical pair norms must be a finite nonnegative square array.');
end
if 3*(32*nnz(graph.pair_norm)+16*numel(ids))>parameters.numerics.max_memory_bytes
    error('cce_partition:memory',...
          'Partition graph and candidate storage exceed numerics.max_memory_bytes.');
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
            'pilot_configs','heldout_configs'};
    if ~all(isfield(p,needed))||~tolerance(p.patch_atol)||~tolerance(p.patch_rtol)||...
       ~integer(p.max_rounds,1)||~integer(p.max_probes,1)||...
       ~iscell(p.patches)||isempty(p.patches)||~isstruct(p.pilot_configs)||...
       isempty(p.pilot_configs)||~isstruct(p.heldout_configs)||isempty(p.heldout_configs)
        error('Adaptive grouping requires finite bounds, patches, pilot, and held-out spin_configs.');
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


