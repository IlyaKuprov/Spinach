% Compile cluster families of physical bath groups. Syntax:
%
%             clusters=cce_clusters(spin_system,graph,parameters)
%
% spin_system is returned by cce_compile; parameters supplies the observed and bath spins.
% graph is cce_graph output. parameters.cce.order counts groups, not sites;
% groups is an optional cell partition of graph.bath_ids (default single
% sites). Non-singleton group construction belongs to cce_partition.
% clusters optionally supplies sorted group-ID rows instead of graph
% generation; cluster_cap is Inf or a per-order nonnegative integer row.
% Connected families grow by graph frontiers, never a bath powerset.
%
% Output is sorted by (order,group IDs); proper_subsets are array indices,
% mobius_power is signed int64, dimension includes the observed-spin system,
% and bath_dimension excludes it. Empty-cluster raw data are separate.
% All nonempty subsets are retained, including disconnected subsets, except
% periodic_segments, which keeps only connected segments of one cycle.
% selected_counts, actual_counts, and cap_overrun report closure costs;
% family_id is the exact canonical family string, not a collision-prone
% numeric hash. config_id=0 denotes the shared unconditioned spin_config.
%
% Parameters:
%
%   spin_system, graph - compiled spin system and physical connectivity.
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

function clusters=cce_clusters(spin_system,graph,parameters)

% Validate the physical group and graph contracts
grumble(spin_system,graph,parameters);

% Resolve the complete physical partition
if isfield(parameters.cce,'groups')&&~isempty(parameters.cce.groups)
    groups=parameters.cce.groups;
else
    groups=num2cell(sort(graph.bath_ids));
end
ng=numel(groups);
max_order=min(parameters.cce.order,ng);

% Disjoint families contain single groups or inter-group pairs
if strcmp(parameters.cce.expansion,'disjoint')
    max_order=min(max_order,1);
elseif strcmp(parameters.cce.expansion,'disjoint_inter')
    max_order=min(max_order,2);
end

% Bound group graph and every subsequent family before growing it
budget=spin_system.tols.cce_mem;
prototype={[]};
metadata=whos('prototype');
family_bytes=metadata.bytes/numel(prototype)+8*(max_order+1);
graph_bytes=16*(nnz(graph.adjacency)+nnz(graph.scores)+2*ng+2);
check_storage(ng);

% Map bath positions to group vertices and sum inter-group scores
[~,sites]=ismember([groups{:}],graph.bath_ids);
vertices=cell(size(groups));
for k=1:ng
    vertices{k}=repmat(k,size(groups{k}));
end
group=sparse(sites,[vertices{:}],1,numel(sites),ng);
adjacency=logical(group.'*sparse(graph.adjacency)*group);
scores=group.'*sparse(graph.scores)*group;
adjacency(1:ng+1:end)=false;
scores(1:ng+1:end)=0;

% Enumerate unique connected candidates by one-vertex frontier growth
explicit=isfield(parameters.cce,'clusters')&&~isempty(parameters.cce.clusters);
if explicit
    check_storage(numel(parameters.cce.clusters));
    family=parameters.cce.clusters(:);
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
    if isfield(parameters.cce,'cluster_cap')&&~isscalar(parameters.cce.cluster_cap)
        count=min(count,parameters.cce.cluster_cap(order));
    elseif isfield(parameters.cce,'cluster_cap')&&isfinite(parameters.cce.cluster_cap)
        count=min(count,parameters.cce.cluster_cap);
    end
    selected=[selected;candidates(1:count)]; %#ok<AGROW>
    selected_counts(order)=count;

end

% Close the selected family unless the expansion uses only periodic segments
family=selected;
closure_order=max_order;
if strcmp(parameters.cce.expansion,'periodic_segments')
    closure_order=1;
end
for order=closure_order:-1:2

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

% Index each sorted group list once for subsequent subset and selection lookups
keys=cellfun(@(ids)string(sprintf('%u,',ids)),family);
selected_keys=cellfun(@(ids)string(sprintf('%u,',ids)),selected);
selected_mask=ismember(keys,selected_keys);
index=dictionary(keys,(1:numel(family)).');

% Bound decimal family signatures, including old/new concatenation copies
physical_digits=1;
if (spin_system.comp.nspins>0)
    physical_digits=1+floor(log10(spin_system.comp.nspins));
end
group_digits=1;
if ng>0
    group_digits=1+floor(log10(ng));
end
characters=(physical_digits+1)*numel(parameters.cce.bath_spins)+2*ng+1+...
           (group_digits+1)*sum(cellfun(@numel,family))+2*numel(family);
if ~isfinite(characters)||6*characters>budget
    error('cce_clusters:memory',...
          'Cluster family signature exceeds sys.tols.cce_mem.');
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
if isfield(parameters.cce,'cluster_cap')&&~all(isinf(parameters.cce.cluster_cap))
    cap=parameters.cce.cluster_cap;
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
             'family_id','','selected_counts',selected_counts,...
             'actual_counts',actual_counts,'cap_overrun',cap_overrun);
metadata=whos('entry');
subset_limit=min(numel(family),2^max_order-1);
bytes=6*characters+3*numel(family)*(metadata.bytes+16*subset_limit+...
                                  8*(4*max_order+3*numel(parameters.cce.bath_spins)));
if ~isfinite(bytes)||bytes>budget
    error('cce_clusters:memory',...
          'Cluster records and subset maps exceed sys.tols.cce_mem.');
end
entry.family_id=family_id;
clusters=repmat(entry,numel(family),1);

% Build canonical physical cluster metadata and subset index maps
[~,observed_pos]=ismember(parameters.observed,1:spin_system.comp.nspins);
observed_dimension=prod(spin_system.comp.mults(observed_pos));
for k=1:numel(family)

    % Map group unions to physical Hilbert dimensions
    sites=sort([groups{family{k}}]);
    [~,positions]=ismember(sites,1:spin_system.comp.nspins);
    clusters(k).id=uint64(k);
    clusters(k).groups=family{k};
    clusters(k).sites=sites;
    clusters(k).order=numel(family{k});
    clusters(k).bath_dimension=prod(spin_system.comp.mults(positions));
    clusters(k).dimension=observed_dimension*clusters(k).bath_dimension;
    if ~isfinite(clusters(k).dimension)||clusters(k).dimension>flintmax
        error('cce_clusters:dimension','Hilbert dimension exceeds exact double integer range.');
    end

    % Find each proper subset in the closed family without scanning other clusters
    subsets=[];
    for order=1:clusters(k).order-1
        rows=nchoosek(family{k},order);
        subkeys=strings(size(rows,1),1);
        for n=1:size(rows,1)
            subkeys(n)=string(sprintf('%u,',rows(n,:)));
        end
        if strcmp(parameters.cce.expansion,'periodic_segments')
            subkeys=subkeys(isKey(index,subkeys));
        end
        children=sort(index(subkeys)).';
        subsets=[subsets children]; %#ok<AGROW>
    end
    clusters(k).proper_subsets=subsets;
    if clusters(k).order>1
        clusters(k).direct_children=children;
    end
    clusters(k).inclusion_reason='closure';
    if selected_mask(k)
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
powers=ones(1,numel(clusters),'int64');
for k=numel(clusters):-1:1

    % Propagate each supercluster contribution only to its actual subsets
    term=powers(k);
    subsets=clusters(k).proper_subsets;
    if (term>0&&any(powers(subsets)<intmin('int64')+term))||...
       (term<0&&any(powers(subsets)>intmax('int64')+term))
        error('cce_clusters:overflow','Mobius coefficient exceeds signed int64 range.');
    end
    powers(subsets)=powers(subsets)-term;
    clusters(k).mobius_power=term;

end

% Bound three family/sorting copies and the sparse group graph
function check_storage(count)

    % Check the workspace size against the explicit budget
    bytes=graph_bytes+3*count*family_bytes;
    if ~isfinite(bytes)||bytes>budget
        error('cce_clusters:memory',...
              'Cluster family workspace exceeds sys.tols.cce_mem.');
    end

end

end

% Validate the immutable metadata and explicit combinatorial decisions
function grumble(spin_system,graph,parameters)
if ~isstruct(parameters)||~isscalar(parameters)||~all(isfield(parameters,{'observed','cce'}))||...
   ~isstruct(parameters.cce)||~isfield(parameters.cce,'bath_spins')
    error('cce_clusters:parameters','resolved parameters are required.');
end
if ~isstruct(spin_system)||~isscalar(spin_system)||~isfield(spin_system,'cce')||...
   ~all(ismember([parameters.observed parameters.cce.bath_spins],1:spin_system.comp.nspins))||...
   ~isempty(intersect(parameters.observed,parameters.cce.bath_spins))
    error('cce_clusters:model','observed and bath spins must be distinct spins of the system.');
end
if ~isstruct(graph)||~isscalar(graph)||~all(isfield(graph,{'bath_ids','adjacency','scores'}))
    error('cce_clusters:graph','cce_graph output is required.');
end
n=numel(graph.bath_ids);
if numel(unique(graph.bath_ids))~=n||~all(ismember(graph.bath_ids,parameters.cce.bath_spins))||...
   ~isequal(size(graph.adjacency),[n n])||~islogical(graph.adjacency)||...
   ~isequal(graph.adjacency,graph.adjacency.')||any(diag(graph.adjacency))||...
   ~isequal(size(graph.scores),[n n])||~isreal(graph.scores)||...
   any(~isfinite(nonzeros(graph.scores)))||any(graph.scores<0,'all')||...
   ~isequal(graph.scores,graph.scores.')
    error('cce_clusters:graph','invalid physical graph ordering, adjacency, or scores.');
end
if ~all(isfield(parameters.cce,{'order','expansion'}))||...
   ~isnumeric(parameters.cce.order)||~isreal(parameters.cce.order)||...
   ~isscalar(parameters.cce.order)||~isfinite(parameters.cce.order)||...
   parameters.cce.order<0||parameters.cce.order~=fix(parameters.cce.order)||...
   ~ismember(parameters.cce.expansion,{'standard','partition','disjoint','disjoint_inter','periodic_segments'})
    error('cce_clusters:parameters','nonnegative order and explicit expansion are required.');
end
if isfield(parameters.cce,'groups')&&~isempty(parameters.cce.groups)
    if ~iscell(parameters.cce.groups)
        error('cce_clusters:groups','groups must be a cell partition of physical bath IDs.');
    end
    for k=1:numel(parameters.cce.groups)
        g=parameters.cce.groups{k};
        if ~isnumeric(g)||~isreal(g)||~isrow(g)||isempty(g)||...
           ~issorted(g)||numel(unique(g))~=numel(g)
            error('cce_clusters:groups','group members must be sorted distinct physical IDs.');
        end
    end
    members=[parameters.cce.groups{:}];
    if numel(unique(members))~=numel(members)||~isequal(sort(members),sort(graph.bath_ids))
        error('cce_clusters:groups','groups must partition the retained bath exactly.');
    end
    ng=numel(parameters.cce.groups);
else
    ng=n;
end
if strcmp(parameters.cce.expansion,'periodic_segments')
    [~,~,components,~]=dmperm(speye(n)+graph.adjacency);
    if n<3||any(sum(graph.adjacency,2)~=2)||numel(components)~=2||...
       parameters.cce.order<1||parameters.cce.order>=n||...
       any(isfield(parameters.cce,{'groups','clusters','cluster_cap'}))
        error('cce_clusters:periodicSegments','Periodic segments require all contiguous subsets of one cycle up to a truncated order.');
    end
end
limit=parameters.cce.order;
if strcmp(parameters.cce.expansion,'disjoint')
    limit=min(limit,1);
end
if strcmp(parameters.cce.expansion,'disjoint_inter')
    limit=min(limit,2);
end
if isfield(parameters.cce,'clusters')&&~isempty(parameters.cce.clusters)
    if ~iscell(parameters.cce.clusters)
        error('cce_clusters:explicit','explicit clusters must be cell rows of group IDs.');
    end
    keys=cell(size(parameters.cce.clusters));
    for k=1:numel(parameters.cce.clusters)
        g=parameters.cce.clusters{k};
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
if isfield(parameters.cce,'cluster_cap')
    cap=parameters.cce.cluster_cap;
    if ~isnumeric(cap)||~isreal(cap)||~isrow(cap)||any(isnan(cap))||any(cap<0)||...
       any(isfinite(cap)&cap~=fix(cap))||...
       (~isscalar(cap)&&numel(cap)<min(parameters.cce.order,ng))||...
       (isscalar(cap)&&isfinite(cap)&&min(parameters.cce.order,ng)>1)
        error('cce_clusters:cap','cluster_cap must be Inf or a sufficiently long integer row.');
    end
end
end


