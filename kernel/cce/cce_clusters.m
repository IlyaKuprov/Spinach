% Compile a downward-closed family of physical bath groups. Syntax:
%
%             clusters=cce_clusters(model,graph,parameters)
%
% model supplies physical_ids, bath_ids, observed_ids, and multiplicities.
% graph is cce_graph output. parameters.order counts groups, not sites;
% groups is an optional cell partition of graph.bath_ids (default single
% sites). Non-singleton group construction belongs to cce_partition.
% clusters optionally supplies sorted group-ID rows instead of graph
% generation; cluster_cap is Inf or a per-order nonnegative integer row.
% Connected families grow by graph frontiers, never a bath powerset.
%
% Output is sorted by (order,group IDs); proper_subsets are array indices,
% mobius_power is signed int64, dimension includes the observed-spin system,
% and bath_dimension excludes it. Empty-cluster raw data are separate.
% All nonempty subsets are retained, including disconnected subsets.
% selected_counts, actual_counts, and cap_overrun report closure costs;
% family_id is the exact canonical family string, not a collision-prone
% numeric hash. config_id=0 denotes the shared unconditioned spin_config.
%
% Parameters:
%
%   model, graph - prepared model and physical connectivity.
%   parameters - order, optional groups (default singleton sites), explicit
%                group-ID clusters, and cluster_cap (default Inf); frozen partition groups
%                are passed in by the driver.
%
% Outputs:
%
%   clusters - ordered records with stable IDs, physical sites, order,
%              groups, dimension, proper_subsets, direct_children, mobius_power,
%              connected_components, config_id, inclusion_reason, weight, family_id,
%              selected_counts, actual_counts, and cap_overrun.
%
% See also cce_graph, cce_partition, cce_reduce.
%
% ilya.kuprov@weizmann.ac.il

function clusters=cce_clusters(model,graph,parameters)

% Validate the physical group and graph contracts
grumble(model,graph,parameters);

% Resolve the complete physical partition
if isfield(parameters,'groups')&&~isempty(parameters.groups)
    groups=parameters.groups;
else
    groups=num2cell(sort(graph.bath_ids));
end
ng=numel(groups);
max_order=min(parameters.order,ng);

% Disjoint families contain single groups or inter-group pairs
if strcmp(parameters.expansion,'disjoint')
    max_order=min(max_order,1);
elseif strcmp(parameters.expansion,'disjoint_inter')
    max_order=min(max_order,2);
end

% Bound group graph and every subsequent family before growing it
budget=Inf;
if isfield(parameters,'numerics')&&isfield(parameters.numerics,'max_memory_bytes')
    budget=parameters.numerics.max_memory_bytes;
end
prototype={[]};
metadata=whos('prototype');
family_bytes=metadata.bytes/numel(prototype)+8*(max_order+1);
check_storage(ng);

% Collapse physical graph edges and scores onto group vertices
adjacency=false(ng);
scores=zeros(ng);
for a=1:ng
    for b=a+1:ng

        % Preserve any physical connection between the composite groups
        [~,ia]=ismember(groups{a},graph.bath_ids);
        [~,ib]=ismember(groups{b},graph.bath_ids);
        adjacency(a,b)=any(graph.adjacency(ia,ib),'all');
        adjacency(b,a)=adjacency(a,b);
        scores(a,b)=sum(graph.scores(ia,ib),'all');
        scores(b,a)=scores(a,b);

    end
end

% Enumerate unique connected candidates by one-vertex frontier growth
explicit=isfield(parameters,'clusters')&&~isempty(parameters.clusters);
if explicit
    check_storage(numel(parameters.clusters));
    family=parameters.clusters(:);
else
    family=cell(0,1);
    frontier=(1:ng).';
    for order=1:max_order

        % Store this order before expanding its graph frontier
        check_storage(numel(family)+size(frontier,1));
        family=[family;mat2cell(frontier,ones(size(frontier,1),1),order)]; %#ok<AGROW>
        next=zeros(0,order+1);
        if order<max_order
            for k=1:size(frontier,1)

                % Add every adjacent vertex, then canonicalise duplicates
                neighbours=find(any(adjacency(frontier(k,:),:),1));
                neighbours=setdiff(neighbours,frontier(k,:));
                for vertex=neighbours
                    check_storage(numel(family)+size(next,1)+1);
                    next(end+1,:)=sort([frontier(k,:) vertex]); %#ok<AGROW>
                end

            end
        end
        frontier=unique(next,'rows');

    end
end

% Select capped candidates by descending score and lexicographic ties
selected=cell(0,1);
selected_counts=zeros(1,max_order);
for order=1:max_order

    % A singleton score is its total incident coupling score
    candidates=family(cellfun(@numel,family)==order);
    if isempty(candidates)
        continue;
    end
    rows=vertcat(candidates{:});
    weights=zeros(size(rows,1),1);
    for k=1:size(rows,1)
        if order==1
            weights(k)=sum(scores(rows(k,:),:),'all');
        else
            weights(k)=sum(scores(rows(k,:),rows(k,:)),'all')/2;
        end
    end
    [~,permutation]=sortrows([-weights rows]);
    candidates=candidates(permutation);

    % Apply selection caps before adding all missing subsets
    count=numel(candidates);
    if isfield(parameters,'cluster_cap')&&~isscalar(parameters.cluster_cap)
        count=min(count,parameters.cluster_cap(order));
    elseif isfield(parameters,'cluster_cap')&&isfinite(parameters.cluster_cap)
        count=min(count,parameters.cluster_cap);
    end
    selected=[selected;candidates(1:count)]; %#ok<AGROW>
    selected_counts(order)=count;

end

% Close the selected family by successive immediate-subset deletion
family=selected;
for order=max_order:-1:2

    % Add immediate subsets and remove duplicate group combinations
    rows=family(cellfun(@numel,family)==order);
    for k=1:numel(rows)
        for omitted=1:order
            check_storage(numel(family)+1);
            family{end+1,1}=rows{k}([1:omitted-1 omitted+1:order]); %#ok<AGROW>
        end
    end
    if ~isempty(family)
        padded=zeros(numel(family),max_order+1);
        for k=1:numel(family)
            padded(k,1:numel(family{k})+1)=[numel(family{k}) family{k}];
        end
        [~,indices]=unique(padded,'rows');
        family=family(indices);
    end

end

% Sort by cardinality and lexical IDs and construct exact family identity
padded=zeros(numel(family),max_order+1);
for k=1:numel(family)
    padded(k,1:numel(family{k})+1)=[numel(family{k}) family{k}];
end
[~,permutation]=sortrows(padded);
family=family(permutation);

% Bound decimal family signatures, including old/new concatenation copies
physical_digits=1;
if ~isempty(model.physical_ids)
    physical_digits=1+floor(log10(max(model.physical_ids)));
end
group_digits=1;
if ng>0
    group_digits=1+floor(log10(ng));
end
characters=(physical_digits+1)*numel(model.bath_ids)+2*ng+1+...
           (group_digits+1)*sum(cellfun(@numel,family))+2*numel(family);
if ~isfinite(characters)||6*characters>budget
    error('cce_clusters:memory',...
          'Cluster family signature exceeds numerics.max_memory_bytes.');
end
family_id='';
for k=1:numel(groups)
    family_id=[family_id '[' sprintf('%u,',groups{k}) ']']; %#ok<AGROW>
end
family_id=[family_id ':'];
for k=1:numel(family)
    family_id=[family_id '[' sprintf('%u,',family{k}) ']']; %#ok<AGROW>
end

% Record actual counts after closure and cap overruns
actual_counts=zeros(1,max_order);
for order=1:max_order
    actual_counts(order)=sum(cellfun(@numel,family)==order);
end
cap_overrun=zeros(1,max_order);
if isfield(parameters,'cluster_cap')&&~all(isinf(parameters.cluster_cap))
    cap=parameters.cluster_cap;
    if isscalar(cap)
        cap=repmat(cap,1,max_order);
    end
    cap_overrun=max(0,actual_counts-cap(1:max_order));
end

% Allocate the fixed cluster schema, including empty-family output
entry=struct('id',uint64(0),'groups',[],'sites',[],'order',0,'dimension',1,...
             'bath_dimension',1,'proper_subsets',[],'direct_children',[],...
             'mobius_power',int64(0),'connected_components',{{}},...
             'config_id',uint64(0),'inclusion_reason','','weight',0,...
             'family_id',family_id,'selected_counts',selected_counts,...
             'actual_counts',actual_counts,'cap_overrun',cap_overrun);
metadata=whos('entry');
bytes=3*numel(family)*(metadata.bytes+16*numel(family)+...
      8*(4*max_order+3*numel(model.bath_ids)));
if ~isfinite(bytes)||bytes>budget
    error('cce_clusters:memory',...
          'Cluster records and subset maps exceed numerics.max_memory_bytes.');
end
clusters=repmat(entry,numel(family),1);

% Build canonical physical cluster metadata and subset index maps
[~,observed_pos]=ismember(model.observed_ids,model.physical_ids);
observed_dimension=prod(model.multiplicities(observed_pos));
for k=1:numel(family)

    % Map group unions to physical Hilbert dimensions
    sites=sort([groups{family{k}}]);
    [~,positions]=ismember(sites,model.physical_ids);
    clusters(k).id=uint64(k);
    clusters(k).groups=family{k};
    clusters(k).sites=sites;
    clusters(k).order=numel(family{k});
    clusters(k).bath_dimension=prod(model.multiplicities(positions));
    clusters(k).dimension=observed_dimension*clusters(k).bath_dimension;
    if ~isfinite(clusters(k).dimension)||clusters(k).dimension>flintmax
        error('cce_clusters:dimension','Hilbert dimension exceeds exact double integer range.');
    end

    % Store every proper nonempty subset, not only immediate children
    subset=false(1,k-1);
    for j=1:k-1
        subset(j)=all(ismember(family{j},family{k}));
    end
    clusters(k).proper_subsets=find(subset);
    clusters(k).direct_children=find(subset&[clusters(1:k-1).order]==clusters(k).order-1);
    clusters(k).inclusion_reason='closure';
    if any(cellfun(@(x)isequal(x,family{k}),selected))
        clusters(k).inclusion_reason='selected';
    end
    clusters(k).weight=sum(scores(family{k},family{k}),'all')/2;

    % Label actual graph components without assuming factorisation
    unseen=family{k};
    components=cell(0,1);
    while ~isempty(unseen)
        component=unseen(1);
        changed=true;
        while changed
            expanded=union(component,intersect(family{k},find(any(adjacency(component,:),1))));
            changed=numel(expanded)>numel(component);
            component=expanded;
        end
        components{end+1,1}=component; %#ok<AGROW>
        unseen=setdiff(unseen,component);
    end
    clusters(k).connected_components=components;

end

% Compute signed inclusion-exclusion powers with checked int64 arithmetic
for k=numel(clusters):-1:1

    % Subtract contributions from all already processed superclusters
    power=int64(1);
    for j=k+1:numel(clusters)
        if ismember(k,clusters(j).proper_subsets)
            term=clusters(j).mobius_power;
            if (term>0&&power<intmin('int64')+term)||...
               (term<0&&power>intmax('int64')+term)
                error('cce_clusters:overflow','Mobius coefficient exceeds signed int64 range.');
            end
            power=power-term;
        end
    end
    clusters(k).mobius_power=power;

end

% Bound three family/sorting copies and the dense group graph
function check_storage(count)

    % Check the workspace size against the explicit budget
    bytes=9*ng^2+3*count*family_bytes;
    if ~isfinite(bytes)||bytes>budget
        error('cce_clusters:memory',...
              'Cluster family workspace exceeds numerics.max_memory_bytes.');
    end

end

end

% Validate the immutable metadata and explicit combinatorial decisions
function grumble(model,graph,parameters)
if ~isstruct(model)||~isscalar(model)||~all(isfield(model,{'physical_ids',...
        'bath_ids','observed_ids','multiplicities'}))||...
   numel(model.physical_ids)~=numel(model.multiplicities)||...
   numel(unique(model.physical_ids))~=numel(model.physical_ids)||...
   ~all(ismember([model.observed_ids model.bath_ids],model.physical_ids))||...
   ~isempty(intersect(model.observed_ids,model.bath_ids))||...
   any(~isfinite(model.multiplicities))||any(model.multiplicities<1)||...
   any(model.multiplicities~=fix(model.multiplicities))
    error('cce_clusters:model','physical IDs and integer multiplicities are required.');
end
if ~isstruct(graph)||~isscalar(graph)||~all(isfield(graph,{'bath_ids','adjacency','scores'}))
    error('cce_clusters:graph','cce_graph output is required.');
end
n=numel(graph.bath_ids);
if numel(unique(graph.bath_ids))~=n||~all(ismember(graph.bath_ids,model.bath_ids))||...
   ~isequal(size(graph.adjacency),[n n])||~islogical(graph.adjacency)||...
   ~isequal(graph.adjacency,graph.adjacency.')||any(diag(graph.adjacency))||...
   ~isequal(size(graph.scores),[n n])||~isreal(graph.scores)||...
   any(~isfinite(graph.scores),'all')||any(graph.scores<0,'all')||...
   ~isequal(graph.scores,graph.scores.')
    error('cce_clusters:graph','invalid physical graph ordering, adjacency, or scores.');
end
if ~isstruct(parameters)||~isscalar(parameters)||...
   ~all(isfield(parameters,{'order','expansion'}))||...
   ~isnumeric(parameters.order)||~isreal(parameters.order)||...
   ~isscalar(parameters.order)||~isfinite(parameters.order)||...
   parameters.order<0||parameters.order~=fix(parameters.order)||...
   ~ismember(parameters.expansion,{'standard','partition','disjoint','disjoint_inter'})
    error('cce_clusters:parameters','nonnegative order and explicit expansion are required.');
end
if isfield(parameters,'groups')&&~isempty(parameters.groups)
    if ~iscell(parameters.groups)
        error('cce_clusters:groups','groups must be a cell partition of physical bath IDs.');
    end
    for k=1:numel(parameters.groups)
        g=parameters.groups{k};
        if ~isnumeric(g)||~isreal(g)||~isrow(g)||isempty(g)||...
           ~issorted(g)||numel(unique(g))~=numel(g)
            error('cce_clusters:groups','group members must be sorted distinct physical IDs.');
        end
    end
    members=[parameters.groups{:}];
    if numel(unique(members))~=numel(members)||~isequal(sort(members),sort(graph.bath_ids))
        error('cce_clusters:groups','groups must partition the retained bath exactly.');
    end
    ng=numel(parameters.groups);
else
    ng=n;
end
limit=parameters.order;
if strcmp(parameters.expansion,'disjoint')
    limit=min(limit,1);
end
if strcmp(parameters.expansion,'disjoint_inter')
    limit=min(limit,2);
end
if isfield(parameters,'clusters')&&~isempty(parameters.clusters)
    if ~iscell(parameters.clusters)
        error('cce_clusters:explicit','explicit clusters must be cell rows of group IDs.');
    end
    keys=cell(size(parameters.clusters));
    for k=1:numel(parameters.clusters)
        g=parameters.clusters{k};
        if ~isnumeric(g)||~isreal(g)||~isrow(g)||isempty(g)||~issorted(g)||...
           numel(unique(g))~=numel(g)||any(g<1)||any(g>ng)||any(g~=fix(g))||...
           numel(g)>limit
            error('cce_clusters:explicit','invalid sorted group-ID cluster or order exceeded.');
        end
        keys{k}=sprintf('%u,',g);
    end
    if numel(unique(keys))~=numel(keys)
        error('cce_clusters:explicit','duplicate explicit clusters are not permitted.');
    end
end
if isfield(parameters,'cluster_cap')
    cap=parameters.cluster_cap;
    if ~isnumeric(cap)||~isreal(cap)||~isrow(cap)||any(isnan(cap))||any(cap<0)||...
       any(isfinite(cap)&cap~=fix(cap))||...
       (~isscalar(cap)&&numel(cap)<min(parameters.order,ng))||...
       (isscalar(cap)&&isfinite(cap)&&min(parameters.order,ng)>1)
        error('cce_clusters:cap','cluster_cap must be Inf or a sufficiently long integer row.');
    end
end
end


