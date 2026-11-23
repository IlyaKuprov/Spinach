% Bind a retained subsystem without constructing a full-bath basis.
% Syntax:
%
%   [spin_system,cache]=cce_subsystem(spin_system,parameters,sites,configuration,cache)
%
% spin_system is cce_compile output and parameters are resolved by
% cce_options. sites is an ordered row of physical bath IDs.
% configuration is a scalar structure; bath_ids and bath_rho specify product
% density matrices in matching order. Omitted densities use the declared
% mixed_product or explicit_product declaration; sampling requires an explicit
% configuration. configuration_id, state_id and weight are optional metadata.
% Explicit exterior_ids/exterior_expectations override no physical states:
% they must agree with the corresponding configuration density matrices.
% cache=[] starts a caller-owned single-entry cache, keyed by the entire spin_system,
% ordered sites and configuration. No process-global or persistent cache exists.
% spin_system.cce carries observed-first IDs, mappings, dimensions, local
% Pauli matrices, local bath densities, and exterior mean fields.
% spin_system.bas.formalism is zeeman-hilb or zeeman-liouv.
% A finite cce.meanfield_radius requires observed and exterior coordinates;
% exterior sites within that Angstrom distance of any observed spin contribute.
% Included sites need no geometry for this state-dependent radius selection.
%
% Parameters:
%
%   spin_system   - spin system returned by cce_compile
%   parameters    - parameters resolved by cce_options
%   sites         - ordered row of retained physical bath IDs
%   configuration - scalar product-state and ensemble-metadata structure
%   cache         - empty or the caller-owned cache returned by this service
%
% Outputs:
%
%   spin_system - the input with the retained ordering, dimensions,
%                 product states, local operators, and exterior spin
%                 expectations in its cce subfield
%   cache       - updated single-entry subsystem cache and hit counters
%
% ilya.kuprov@weizmann.ac.il

function [spin_system,cache]=cce_subsystem(spin_system,parameters,sites,configuration,cache)

% Reject malformed ownership and configuration before consulting the cache
grumble(spin_system,parameters,sites,configuration,cache);

% An exact key avoids hash collisions and stale spin_system or representation data
signature=struct('spin_system',spin_system,'parameters',parameters,'sites',sites,'configuration',configuration);
if ~isempty(cache)&&isequaln(cache.signature,signature)
    spin_system=cache.spin_system;
    cache.hits=cache.hits+1;
    return
end

% Preserve observed-first tensor-product order and explicit physical mappings
ids=[parameters.observed sites];
dims=spin_system.comp.mults(ids);
spin_system.cce.sites=sites;
spin_system.cce.ids=ids;
spin_system.cce.physical_to_local=zeros(1,spin_system.comp.nspins);
spin_system.cce.physical_to_local(ids)=1:numel(ids);
spin_system.cce.dimension=prod(dims);
spin_system.cce.observed_dimension=prod(spin_system.comp.mults(parameters.observed));
spin_system.cce.bath_dimension=prod(spin_system.comp.mults(sites));
spin_system.cce.groups={};
if isfield(parameters.cce,'groups')
    spin_system.cce.groups=cellfun(@(group)group(ismember(group,ids)),...
                                  parameters.cce.groups,'UniformOutput',false);
end
spin_system.cce.local_operators=cellfun(@pauli,num2cell(dims),'UniformOutput',false);

% Bind a declared product bath density, never an inferred population
bath_ids=parameters.cce.bath_spins;
if isfield(configuration,'bath_rho')
    bath_ids=configuration.bath_ids;
    densities=configuration.bath_rho;
elseif isfield(parameters,'bath_state')
    densities=parameters.bath_state;
else
    densities=arrayfun(@(id)speye(spin_system.comp.mults(id))/spin_system.comp.mults(id),...
                       bath_ids,'UniformOutput',false);
end
expectations=zeros(numel(bath_ids),3);
for n=1:numel(bath_ids)

    % Compute physical Cartesian expectations from each validated local density
    d=spin_system.comp.mults(bath_ids(n));
    rho=densities{n};
    operators=pauli(d);
    expectations(n,:)=[real(trace(rho*operators.x)) real(trace(rho*operators.y))...
                      real(trace(rho*operators.z))];

end
if isfield(configuration,'expectations')&&...
   ~close_values(configuration.expectations,expectations,spin_system.tols)
    error('cce_subsystem:expectations','Spin expectations disagree with the supplied density.');
end

% Use the same exterior assignment for every retained cluster
[~,local]=ismember(sites,bath_ids);
spin_system.cce.bath_rho=densities(local);
exterior=~ismember(bath_ids,sites);
exterior_ids=bath_ids(exterior);
exterior_values=expectations(exterior,:);
if isfield(configuration,'exterior_ids')
    supplied=configuration.exterior_ids;
    [present,indices]=ismember(supplied,exterior_ids);
    if ~all(present)||numel(indices)~=numel(exterior_ids)||...
       ~close_values(configuration.exterior_expectations,...
                     exterior_values(indices,:),spin_system.tols)
        error('cce_subsystem:exterior','Exterior fields must match the common bath configuration.');
    end
end
spin_system.cce.configuration=configuration;
spin_system.cce.configuration.bath_ids=bath_ids;
spin_system.cce.configuration.bath_rho=densities;
spin_system.cce.configuration.expectations=expectations;
spin_system.cce.configuration.exterior_ids=[];
spin_system.cce.configuration.exterior_expectations=zeros(0,3);
if strcmp(parameters.cce.meanfield,'static')
    if isfield(parameters.cce,'meanfield_radius')&&...
       isfinite(parameters.cce.meanfield_radius)&&~isempty(exterior_ids)
        distance=inf(1,numel(exterior_ids));
        for id=parameters.observed
            delta=vertcat(spin_system.inter.coordinates{exterior_ids})-spin_system.inter.coordinates{id};
            distance=min(distance,sqrt(sum(delta.^2,2)).');
        end
        keep=distance<=parameters.cce.meanfield_radius;
        exterior_ids=exterior_ids(keep);
        exterior_values=exterior_values(keep,:);
    end
    spin_system.cce.configuration.exterior_ids=exterior_ids;
    spin_system.cce.configuration.exterior_expectations=exterior_values;
end

% Cache only this explicit immutable view; generator caches belong downstream
hits=0;
misses=0;
if ~isempty(cache)
    hits=cache.hits;
    misses=cache.misses;
end
cache=struct('signature',signature,'spin_system',spin_system,...
             'hits',hits,'misses',misses+1);

end

% Compare finite physical expectations with the declared numerical tolerance
function answer=close_values(left,right,tols)
answer=isnumeric(left)&&isequal(size(left),size(right))&&...
       all(isfinite(left),'all')&&...
       norm(left-right,'fro')<=tols.cce_atol+tols.cce_rtol*norm(right,'fro');
end

% Validate the complete public service boundary
function grumble(spin_system,parameters,sites,configuration,cache)
if ~isstruct(spin_system)||~isscalar(spin_system)||~isfield(spin_system,'cce')
    error('cce_subsystem:model','model must be returned by cce_compile.');
end
if ~isnumeric(sites)||~isreal(sites)||(~isempty(sites)&&~isrow(sites))||...
   any(~ismember(sites,parameters.cce.bath_spins))||numel(unique(sites))~=numel(sites)
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
       ~isequal(sort(ids),sort(parameters.cce.bath_spins))
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
   ~all(isfield(cache,{'signature','spin_system','hits','misses'})))
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
        d=spin_system.comp.mults(bath_ids(n));
        rho=densities{n};
        tolerance=spin_system.tols.cce_atol+spin_system.tols.cce_rtol*norm(rho,'fro');
        if ~isnumeric(rho)||~isequal(size(rho),[d d])||any(~isfinite(rho),'all')||...
           norm(rho-rho','fro')>tolerance||abs(trace(rho)-1)>tolerance||...
           min(real(eig(full((rho+rho')/2))))<-tolerance
            error('cce_subsystem:density','Bath densities must be positive, Hermitian and trace one.');
        end
    end
elseif isfield(parameters,'bath_samples')
    error('cce_subsystem:configuration','Sampling requires explicit bath density matrices.');
end
if isfield(configuration,'exterior_ids')
    supplied=configuration.exterior_ids;
    if ~isnumeric(supplied)||~isreal(supplied)||...
       (~isempty(supplied)&&~isrow(supplied))||numel(unique(supplied))~=numel(supplied)
        error('cce_subsystem:exterior','Exterior IDs must be a distinct physical-ID row.');
    end
end
if strcmp(parameters.cce.meanfield,'static')&&...
   isfield(parameters.cce,'meanfield_radius')&&...
   isfinite(parameters.cce.meanfield_radius)
    exterior=setdiff(parameters.cce.bath_spins,sites);
    if ~isempty(exterior)&&...
       any(cellfun(@isempty,spin_system.inter.coordinates([parameters.observed exterior])))
        error('cce:coordinates','cce.meanfield_radius requires finite observed and exterior coordinates.');
    end
end
end


