% Bind a retained subsystem without constructing a full-bath basis.
% Syntax:
%
%   [cluster_system,cache]=cce_subsystem(model,sites,configuration,cache)
%
% model is cce_compile output. sites is an ordered row of physical bath IDs.
% configuration is a scalar structure; bath_ids and bath_rho specify product
% density matrices in matching order. Omitted densities use the model's
% mixed_product or explicit_product declaration; sampling requires an explicit
% configuration. configuration_id, state_id and weight are optional metadata.
% Explicit exterior_ids/exterior_expectations override no physical states:
% they must agree with the corresponding configuration density matrices.
% cache=[] starts a caller-owned single-entry cache, keyed by the entire model,
% ordered sites and configuration. No process-global or persistent cache exists.
% cluster_system carries observed-first IDs, mappings, multiplicities, isotope
% labels, local Pauli matrices, local bath densities, and exterior mean fields.
% formalism is zeeman-hilb or zeeman-liouv; Hilbert dimensions are always explicit.
% A finite cce.meanfield_radius requires observed and exterior coordinates;
% exterior sites within that Angstrom distance of any observed spin contribute.
% Included sites need no geometry for this state-dependent radius selection.
%
% Parameters:
%
%   model         - canonical physical model returned by cce_compile
%   sites         - ordered row of retained physical bath IDs
%   configuration - scalar product-state and ensemble-metadata structure
%   cache         - empty or the caller-owned cache returned by this service
%
% Outputs:
%
%   cluster_system - retained tensor ordering, dimensions, product states,
%                    local operators, and exterior spin expectations
%   cache          - updated single-entry subsystem cache and hit counters
%
% ilya.kuprov@weizmann.ac.il

function [cluster_system,cache]=cce_subsystem(model,sites,configuration,cache)

% Reject malformed ownership and configuration before consulting the cache
grumble(model,sites,configuration,cache);

% An exact key avoids hash collisions and stale model or representation data
signature=struct('model',model,'sites',sites,'configuration',configuration);
if ~isempty(cache)&&isequaln(cache.signature,signature)
    cluster_system=cache.cluster_system;
    cache.hits=cache.hits+1;
    return
end

% Preserve observed-first tensor-product order and explicit physical mappings
ids=[model.observed_ids sites];
dims=model.multiplicities(ids);
cluster_system.kind='cce_subsystem';
cluster_system.model=model;
cluster_system.formalism=model.formalism;
cluster_system.sites=sites;
cluster_system.physical_ids=ids;
cluster_system.local_to_physical=ids;
cluster_system.physical_to_local=zeros(size(model.physical_ids));
cluster_system.physical_to_local(ids)=1:numel(ids);
cluster_system.multiplicities=dims;
cluster_system.dimension=prod(dims);
cluster_system.observed_dimension=prod(model.multiplicities(model.observed_ids));
cluster_system.bath_dimension=prod(model.multiplicities(sites));
cluster_system.isotopes=model.isotopes(ids);
cluster_system.coordinates_A=model.coordinates_A(ids,:);
cluster_system.assumptions=model.assumptions;
cluster_system.groups=cellfun(@(group)group(ismember(group,ids)),...
                              model.groups,'UniformOutput',false);
cluster_system.local_operators=cellfun(@pauli,num2cell(dims),'UniformOutput',false);

% Bind a declared product bath density, never an inferred population
bath_ids=model.bath_ids;
state=model.parameters.cce.bath_state;
if isfield(configuration,'bath_rho')
    bath_ids=configuration.bath_ids;
    densities=configuration.bath_rho;
elseif strcmp(state.mode,'mixed_product')
    densities=arrayfun(@(id)speye(model.multiplicities(id))/model.multiplicities(id),...
                       bath_ids,'UniformOutput',false);
elseif strcmp(state.mode,'explicit_product')
    densities=state.rho;
else
    error('cce_subsystem:configuration','Sampling requires explicit bath density matrices.');
end
expectations=zeros(numel(bath_ids),3);
for n=1:numel(bath_ids)

    % Compute physical Cartesian expectations from each validated local density
    d=model.multiplicities(bath_ids(n));
    rho=densities{n};
    operators=pauli(d);
    expectations(n,:)=[real(trace(rho*operators.x)) real(trace(rho*operators.y))...
                      real(trace(rho*operators.z))];

end
if isfield(configuration,'expectations')&&...
   ~close_values(configuration.expectations,expectations,model.numerics)
    error('cce_subsystem:expectations','Spin expectations disagree with the supplied density.');
end

% Use the same exterior assignment for every retained cluster
[~,local]=ismember(sites,bath_ids);
cluster_system.bath_rho=densities(local);
exterior=~ismember(bath_ids,sites);
exterior_ids=bath_ids(exterior);
exterior_values=expectations(exterior,:);
if isfield(configuration,'exterior_ids')
    supplied=configuration.exterior_ids;
    [present,indices]=ismember(supplied,exterior_ids);
    if ~all(present)||numel(indices)~=numel(exterior_ids)||...
       ~close_values(configuration.exterior_expectations,...
                     exterior_values(indices,:),model.numerics)
        error('cce_subsystem:exterior','Exterior fields must match the common bath configuration.');
    end
end
cluster_system.configuration=configuration;
cluster_system.configuration.bath_ids=bath_ids;
cluster_system.configuration.bath_rho=densities;
cluster_system.configuration.expectations=expectations;
cluster_system.configuration.exterior_ids=[];
cluster_system.configuration.exterior_expectations=zeros(0,3);
if strcmp(model.parameters.cce.meanfield,'static')
    if isfield(model.parameters.cce,'meanfield_radius')&&...
       isfinite(model.parameters.cce.meanfield_radius)&&~isempty(exterior_ids)
        coordinates=model.coordinates_A([model.observed_ids exterior_ids],:);
        if any(~isfinite(coordinates),'all')
            error('cce:coordinates','cce.meanfield_radius requires finite observed and exterior coordinates.');
        end
        distance=inf(1,numel(exterior_ids));
        for id=model.observed_ids
            delta=model.coordinates_A(exterior_ids,:)-model.coordinates_A(id,:);
            distance=min(distance,sqrt(sum(delta.^2,2)).');
        end
        keep=distance<=model.parameters.cce.meanfield_radius;
        exterior_ids=exterior_ids(keep);
        exterior_values=exterior_values(keep,:);
    end
    cluster_system.configuration.exterior_ids=exterior_ids;
    cluster_system.configuration.exterior_expectations=exterior_values;
end

% Cache only this explicit immutable view; generator caches belong downstream
hits=0;
misses=0;
if ~isempty(cache)
    hits=cache.hits;
    misses=cache.misses;
end
cache=struct('signature',signature,'cluster_system',cluster_system,...
             'hits',hits,'misses',misses+1);

end

% Compare finite physical expectations with the declared numerical tolerance
function answer=close_values(left,right,numerics)
answer=isnumeric(left)&&isequal(size(left),size(right))&&...
       all(isfinite(left),'all')&&...
       norm(left-right,'fro')<=numerics.atol+numerics.rtol*norm(right,'fro');
end

% Validate the complete public service boundary
function grumble(model,sites,configuration,cache)
if ~isstruct(model)||~isscalar(model)||~isfield(model,'kind')||...
   ~strcmp(model.kind,'cce_model')
    error('cce_subsystem:model','model must be returned by cce_compile.');
end
if ~isnumeric(sites)||~isreal(sites)||(~isempty(sites)&&~isrow(sites))||...
   any(~ismember(sites,model.bath_ids))||numel(unique(sites))~=numel(sites)
    error('cce_subsystem:sites','sites must be distinct physical bath IDs in the requested order.');
end
if ~isstruct(configuration)||~isscalar(configuration)
    error('cce_subsystem:configuration','configuration must be a scalar structure.');
end
allowed={'bath_ids','bath_rho','expectations','exterior_ids','exterior_expectations',...
         'configuration_id','state_id','weight','backend'};
if ~isempty(setdiff(fieldnames(configuration),allowed))
    error('cce_subsystem:configuration','Unknown configuration field.');
end
if xor(isfield(configuration,'bath_ids'),isfield(configuration,'bath_rho'))
    error('cce_subsystem:configuration','bath_ids and bath_rho must be supplied together.');
end
if isfield(configuration,'bath_ids')
    ids=configuration.bath_ids;
    if ~isnumeric(ids)||(~isempty(ids)&&~isrow(ids))||...
       ~isequal(sort(ids),sort(model.bath_ids))
        error('cce_subsystem:configuration','bath_ids must list the selected bath exactly once.');
    end
end
if xor(isfield(configuration,'exterior_ids'),isfield(configuration,'exterior_expectations'))
    error('cce_subsystem:configuration','Exterior IDs and expectations must be supplied together.');
end
if isfield(configuration,'backend')&&...
   ~ismember(configuration.backend,{'auto','direct','stock'})
    error('cce_subsystem:backend','Unknown Hamiltonian backend.');
end
if ~isempty(cache)&&(~isstruct(cache)||~isscalar(cache)||...
   ~all(isfield(cache,{'signature','cluster_system','hits','misses'})))
    error('cce_subsystem:cache','cache must be empty or returned by this service.');
end
if isfield(configuration,'bath_rho')
    bath_ids=configuration.bath_ids;
    densities=configuration.bath_rho;
    if ~iscell(densities)||(~isempty(densities)&&~isrow(densities))||...
       numel(densities)~=numel(bath_ids)
        error('cce_subsystem:density','Supply one density matrix per bath ID.');
    end
    for n=1:numel(bath_ids)
        d=model.multiplicities(bath_ids(n));
        rho=densities{n};
        tolerance=model.numerics.atol+model.numerics.rtol*norm(rho,'fro');
        if ~isnumeric(rho)||~isequal(size(rho),[d d])||any(~isfinite(rho),'all')||...
           norm(rho-rho','fro')>tolerance||abs(trace(rho)-1)>tolerance||...
           min(real(eig(full((rho+rho')/2))))<-tolerance
            error('cce_subsystem:density','Bath densities must be positive, Hermitian and trace one.');
        end
    end
elseif ismember(model.parameters.cce.bath_state.mode,{'sampled','interlaced'})
    error('cce_subsystem:configuration','Sampling requires explicit bath density matrices.');
end
if isfield(configuration,'exterior_ids')
    supplied=configuration.exterior_ids;
    if ~isnumeric(supplied)||~isreal(supplied)||...
       (~isempty(supplied)&&~isrow(supplied))||numel(unique(supplied))~=numel(supplied)
        error('cce_subsystem:exterior','Exterior IDs must be a distinct physical-ID row.');
    end
end
end


