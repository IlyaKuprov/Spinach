% Bath connectivity for the cluster expansion. Syntax:
%
%                  graph=cce_graph(spin_system,parameters)
%
% Bath spins are joined when create kept a coupling between them, so the
% connectivity follows sys.tols.prox_cutoff and sys.tols.inter_cutoff. In
% generalised and second-order projected expansions, bath spins coupled to
% the observed spins are also joined to one another, because the observed
% spins mediate their correlation; parameters.cce.mediated=false drops
% those edges. Noise calculations join the spins that support the noise
% operators. parameters.cce.bath_radius (Angstrom, default Inf) restricts
% the bath to spins within that distance of an observed spin.
%
% Parameters:
%
%   spin_system - spin system returned by cce_compile
%
%   parameters  - observed and cce.bath_spins physical spin IDs, with
%                 cce.approximation, cce.second_order, cce.mediated,
%                 and optional cce.bath_radius
%
% Outputs:
%
%   graph - structure with bath_ids, sparse logical adjacency, pair_norm
%           (two-norms of the pair coupling operators, rad/s), scores
%           (equal to pair_norm, used to rank clusters), observed_mediated,
%           and omitted_bath_ids; matrix ordering follows bath_ids
%
% See also cce_clusters, cce_partition.
%
% ilya.kuprov@weizmann.ac.il

function graph=cce_graph(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Bath spin IDs in the caller's order
ids=parameters.cce.bath_spins;
n=numel(ids);

% Two-norms of the pair coupling operators
ends=reshape([spin_system.cce.pairs.ids],2,[]);
weights=zeros(1,numel(spin_system.cce.pairs));
for k=1:numel(spin_system.cce.pairs)

    % Stock spin matrices at the multiplicities of the two spins
    A=pauli(spin_system.comp.mults(ends(1,k)));
    B=pauli(spin_system.comp.mults(ends(2,k)));
    A={A.x,A.y,A.z};
    B={B.x,B.y,B.z};

    % Bilinear coupling operator
    T=spin_system.cce.pairs(k).tensor;
    C=0;
    for u=1:3
        for v=1:3
            C=C+T(u,v)*kron(A{u},B{v});
        end
    end
    weights(k)=norm(full(C),2);

end

% Coupling norms between bath spins
[~,bath]=ismember(ends,ids);
internal=all(bath>0,1);
a=bath(1,internal);
b=bath(2,internal);
pair_norm=sparse([a b],[b a],[weights(internal) weights(internal)],n,n);

% Bath spins coupled to the observed spins
observed=any(ismember(ends,parameters.observed),1);
attachments=sum(bath(:,observed&(weights>0)),1);
attached=sparse(attachments(attachments>0),1,true,n,1);
mediated=logical(attached*attached.');
mediated(1:n+1:end)=false;

% Join the bath spins that create left coupled
adjacency=pair_norm>0;

% Keep the correlations mediated by the observed spins
second_order=strcmp(parameters.cce.approximation,'projected')&&...
             parameters.cce.second_order.enabled;
if ((strcmp(parameters.cce.approximation,'generalized')&&...
     ~isfield(parameters,'noise'))||second_order)&&parameters.cce.mediated
    adjacency=adjacency|mediated;
end

% Retain disconnected moments and finite-sample cross-site covariance
if isfield(parameters,'noise')&&...
   (~parameters.noise.connected||isfield(parameters,'bath_samples'))
    if isfield(parameters.noise,'terms')
        support=[parameters.noise.terms.support];
    else
        support=[parameters.noise.descriptors.spins];
    end
    support=ismember(ids,support);
    adjacency(support,support)=true;
    adjacency(1:n+1:end)=false;
end

% Restrict the bath by distance from the nearest observed spin
keep=true(1,n);
if isfield(parameters.cce,'bath_radius')&&isfinite(parameters.cce.bath_radius)
    distance=inf(n,1);
    for k=parameters.observed
        distance=min(distance,vecnorm(vertcat(spin_system.inter.coordinates{ids})-...
                                      spin_system.inter.coordinates{k},2,2));
    end
    keep=distance.'<=parameters.cce.bath_radius;
end

% Return the graph in the physical bath ordering
graph=struct('bath_ids',ids(keep),'adjacency',sparse(adjacency(keep,keep)),...
             'scores',pair_norm(keep,keep),'pair_norm',pair_norm(keep,keep),...
             'observed_mediated',sparse(mediated(keep,keep)),...
             'omitted_bath_ids',ids(~keep));

end

% Consistency enforcement
function grumble(spin_system,parameters)
if ~isstruct(parameters)||~isscalar(parameters)||~all(isfield(parameters,{'observed','cce'}))||...
   ~all(isfield(parameters.cce,{'bath_spins','approximation','mediated'}))
    error('cce_graph:parameters','resolved cluster expansion parameters are required.');
end
if ~islogical(parameters.cce.mediated)||~isscalar(parameters.cce.mediated)
    error('cce_graph:mediated','parameters.cce.mediated must be a logical scalar.');
end
if ~isstruct(spin_system)||~isscalar(spin_system)||~isfield(spin_system,'cce')||...
   ~isfield(spin_system.cce,'pairs')
    error('cce_graph:model','incomplete scalar spin_system structure.');
end
ids=1:spin_system.comp.nspins;
for name={parameters.observed,parameters.cce.bath_spins}
    values=name{1};
    if ~isnumeric(values)||~isreal(values)||(~isrow(values)&&~isempty(values))||...
       numel(unique(values))~=numel(values)||~all(ismember(values,ids))
        error('cce_graph:ids','observed and bath IDs must be unique physical ID rows.');
    end
end
if ~isempty(intersect(parameters.observed,parameters.cce.bath_spins))||isempty(parameters.observed)
    error('cce_graph:ids','observed IDs must be nonempty and disjoint from bath IDs.');
end
if ~isstruct(spin_system.cce.pairs)||(~isempty(spin_system.cce.pairs)&&~all(isfield(spin_system.cce.pairs,{'ids','tensor'})))
    error('cce_graph:pairs','pairs must be records ids,tensor.');
end
pairs=zeros(numel(spin_system.cce.pairs),2);
for k=1:numel(spin_system.cce.pairs)
    p=spin_system.cce.pairs(k);
    if ~isnumeric(p.ids)||~isrow(p.ids)||numel(p.ids)~=2||...
       ~all(ismember(p.ids,ids))||p.ids(1)==p.ids(2)||...
       ~isnumeric(p.tensor)||~isreal(p.tensor)||~isequal(size(p.tensor),[3 3])||...
       any(~isfinite(p.tensor),'all')
        error('cce_graph:pairs','invalid physical pair or Cartesian rad/s tensor.');
    end
    pairs(k,:)=sort(p.ids);
end
if size(unique(pairs,'rows'),1)~=size(pairs,1)
    error('cce_graph:pairs','unordered pairs must already be combined.');
end
if isfield(parameters.cce,'bath_radius')&&(~isnumeric(parameters.cce.bath_radius)||...
   ~isreal(parameters.cce.bath_radius)||~isscalar(parameters.cce.bath_radius)||...
   isnan(parameters.cce.bath_radius)||parameters.cce.bath_radius<0)
    error('cce_graph:radius','bath_radius must be nonnegative.');
end
if isfield(parameters.cce,'bath_radius')&&isfinite(parameters.cce.bath_radius)&&...
   any(cellfun(@isempty,spin_system.inter.coordinates))
    error('cce_graph:coordinates','bath_radius requires coordinates for every spin.');
end
end


