% Projected product-state averages for interlaced CCE (W12). Syntax:
%
%       raw=cce_interlace(spin_system,cluster,clusters,ensemble,pulse_sequence,parameters)
%
% parameters.cce.interlaced contains max_states (explicit enumeration budget),
% reference_state (ensemble spin_config index), and optionally condition_sites
% (cell array of physical-ID sets). Otherwise K is the union of cylinders
% associated with included superclusters of cluster. Denominator subsets
% MUST reuse raw.condition_sites, not their own default cylinder family.
% The reference bath densities must agree with their basis assignment.
% The returned raw value is the probability-weighted conditional average;
% configurations are unique and weights are conditional probabilities.
% The distribution signature records assignments and weights exactly.
%
% Parameters:
%
%   spin_system, cluster, clusters, ensemble, pulse_sequence, parameters -
%   compiled spin system, selected cluster, cluster family, complete sampled
%   ensemble, experiment function handle, and interlaced settings.
%   parameters.cce.interlaced - max_states, reference_state (ensemble spin_config
%   index), optional condition_sites cell array and raw_only flag.
%
% Outputs:
%
%   raw - conditional raw cluster observable with condition_sites and
%         conditional-average diagnostics for consistent denominator reduction.
%
% Native controls are nested in parameters.cce. Native raw values stay
% conditional averages; factor fields hold the separately reduced quotient.
% Assignment probabilities use shifted logarithmic products, and numerical
% residuals retain the worst constituent solve. raw_only leaves factors empty.
%
% See also cce_reduce.
% Reference: https://doi.org/10.1103/PhysRevB.86.035452
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_interlace(spin_system,cluster,clusters,ensemble,pulse_sequence,parameters)

% Validate the native conditional ensemble and explicit experiment callback
grumble(spin_system,cluster,clusters,ensemble,pulse_sequence,parameters);

% Reuse the same callback and conditioning event in every denominator
raw=native_interlace(spin_system,cluster,clusters,ensemble,pulse_sequence,parameters);

end

% Form the conditional average and its separately retained irreducible factor
function raw=native_interlace(spin_system,cluster,clusters,ensemble,pulse_sequence,parameters)

% Select the validated state budget and common reference assignment
spec=parameters.cce.interlaced;
reference=ensemble.spin_config(spec.reference_state);

% Resolve the supplied cylinders or the selected supercluster family
if isfield(spec,'condition_sites')
    cylinders=spec.condition_sites;
else
    cylinders={};
    for k=1:numel(clusters)
        parent=clusters(k);
        if all(ismember(cluster.sites,parent.sites))
            cylinders{end+1}=parent.sites; %#ok<AGROW>
        end
    end
    if isempty(cylinders)
        cylinders={cluster.sites};
    end
end

% Remove contained cylinders before charging the explicit enumeration budget
keep=true(size(cylinders));
for j=1:numel(cylinders)
    for k=1:numel(cylinders)
        if j~=k&&all(ismember(cylinders{j},cylinders{k}))&&...
           (numel(cylinders{j})<numel(cylinders{k})||j>k)
            keep(j)=false;
        end
    end
end
cylinders=cylinders(keep);
nb=numel(parameters.cce.bath_spins);
dimensions=spin_system.comp.mults(parameters.cce.bath_spins);
counts=cellfun(@(sites)prod(dimensions(ismember(parameters.cce.bath_spins,sites))),cylinders);
count=sum(counts);
if count>spec.max_states
    error('cce_interlace:budget','Cylinder enumeration exceeds interlaced.max_states.');
end
bytes=8*count*(4*nb+5);
if ~isfinite(bytes)||bytes>spin_system.tols.cce_mem
    error('cce_interlace:memory','Conditioning assignments and sorting workspace exceed the memory budget.');
end

% Enumerate unique product-state assignments inside the cylinder union
assignments=zeros(count,nb);
position=0;
for k=1:numel(cylinders)
    varied=find(ismember(parameters.cce.bath_spins,cylinders{k}));
    for n=1:counts(k)
        position=position+1;
        assignments(position,:)=reference.assignment;
        index=n-1;
        for j=fliplr(varied)
            assignments(position,j)=mod(index,dimensions(j))+1;
            index=floor(index/dimensions(j));
        end
    end
end
assignments=unique(assignments,'rows','sorted');
count=size(assignments,1);

% Form conditional weights with a common logarithmic scale
logweights=zeros(count,1);
for j=1:nb
    logweights=logweights+log(ensemble.probabilities{j}(assignments(:,j)));
end
shift=max(logweights);
if ~isfinite(shift)
    error('cce_interlace:zero_mass','Conditioning event has zero probability.');
end
weights=exp(logweights-shift);
logmass=shift+log(sum(weights));
weights=weights/sum(weights);
mass=exp(logmass);
retained=find(weights>0);

% Average raw records before forming any irreducible quotient
for n=reshape(retained,1,[])
    configuration=reference;
    configuration.assignment=assignments(n,:);
    for j=1:nb
        vector=ensemble.bases{j}(:,assignments(n,j));
        configuration.bath_rho{j}=vector*vector';
    end
    item=cce_evaluate(spin_system,cluster,configuration,pulse_sequence,parameters);
    if n==retained(1)
        raw=item;
        raw.value(:)=0;
        raw.abs_error_bound(:)=0;
        raw.valid(:)=true;
    end
    raw.value=raw.value+weights(n)*item.value;
    raw.abs_error_bound=raw.abs_error_bound+weights(n)*item.abs_error_bound;
    raw.valid=raw.valid&item.valid;
    for field={'eigen_residual','orthogonality_residual','unitarity_residual'}
        if isfield(item,field{1})
            raw.(field{1})=max(raw.(field{1}),item.(field{1}));
        end
    end
end
raw.condition_sites=cylinders;
raw.assignment_count=numel(retained);
raw.assignments=assignments;
raw.weights=weights;
raw.condition_mass=mass;
raw.log_condition_mass=logmass;
raw.distribution_signature=struct('assignments',assignments,'weights',weights);
raw.factor=[];
raw.factor_valid=[];
raw.factor_abs_error_bound=[];
raw.reference=raw.value;
raw.reference_error=raw.abs_error_bound;
raw.full_bath=isequal(sort(cluster.sites),sort(parameters.cce.bath_spins));
if isfield(spec,'raw_only')&&spec.raw_only
    return;
end

% Every denominator receives exactly the parent's conditioning event
if ~isempty(cluster.sites)&&~raw.full_bath
    local=parameters;
    local.cce.interlaced.condition_sites=cylinders;
    local.cce.interlaced.raw_only=true;
    family=clusters;
    included=arrayfun(@(member)all(ismember(member.sites,cluster.sites)),family);
    family=family(included);
    records=cell(1,numel(family)+1);
    empty=struct('id',0,'sites',[]);
    records{1}=native_interlace(spin_system,empty,clusters,ensemble,pulse_sequence,local);
    for k=1:numel(family)
        if isequal(family(k).sites,cluster.sites), records{k+1}=raw;
        else
            records{k+1}=native_interlace(spin_system,family(k),clusters,ensemble,pulse_sequence,local);
        end
    end
    factors=records{1}.value;
    errors=records{1}.abs_error_bound;
    for k=1:numel(family)
        item=records{k+1};
        denominator=records{1}.value;
        denominator_error=records{1}.abs_error_bound;
        valid=records{1}.valid&item.valid;
        for j=family(k).proper_subsets
            position=find([family.id]==clusters(j).id,1);
            if isempty(position)
                continue;
            end
            factor=records{position+1};
            denominator_error=abs(denominator).*factor.abs_error_bound+...
                              abs(factor.value).*denominator_error+...
                              denominator_error.*factor.abs_error_bound;
            denominator=denominator.*factor.value;
            valid=valid&factor.valid;
        end
        resolved=valid&abs(denominator)>denominator_error+spin_system.tols.cce_atol+...
                 spin_system.tols.cce_rtol*abs(denominator);
        value=NaN(size(item.value));
        uncertainty=Inf(size(item.value));
        value(resolved)=item.value(resolved)./denominator(resolved);
        uncertainty(resolved)=(item.abs_error_bound(resolved)+...
            abs(value(resolved)).*denominator_error(resolved))./...
            (abs(denominator(resolved))-denominator_error(resolved));
        item.value=value;
        item.abs_error_bound=uncertainty;
        item.valid=resolved&isfinite(value);
        records{k+1}=item;
    end
    position=find([family.id]==cluster.id,1);
    factor=records{position+1};
    raw.factor=factor.value;
    raw.factor_valid=factor.valid;
    raw.factor_abs_error_bound=factor.abs_error_bound;
    raw.reference=factors;
    raw.reference_error=errors;
end
raw.kind='irreducible';
if raw.full_bath
    raw.kind='raw_full_bath';
end
end

% Validate the resolved ensemble contract without a method-name fallback
function grumble(spin_system,cluster,clusters,ensemble,pulse_sequence,parameters)
if ~isstruct(spin_system)||~isfield(spin_system,'cce')||...
   ~isstruct(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   (~isempty(clusters)&&(~isstruct(clusters)||~all(isfield(clusters,{'id','sites'}))))||...
   ~isstruct(ensemble)||~all(isfield(ensemble,{'spin_config','bases','probabilities'}))||...
   ~isa(pulse_sequence,'function_handle')||...
   ~isstruct(parameters)||~isfield(parameters,'cce')||...
   ~isfield(parameters.cce,'interlaced')
    error('cce_interlace:input','Native spin_system, cluster, ensemble and experiment handle required.');
end
if isfield(parameters.cce,'checkpoint')
    error('cce_options:checkpoint','Interlaced evaluations cannot own a checkpoint ledger.');
end
spec=parameters.cce.interlaced;
if ~isstruct(spec)||~all(isfield(spec,{'reference_state','max_states'}))||...
   ~isscalar(spec.reference_state)||~isreal(spec.reference_state)||...
   ~isfinite(spec.reference_state)||spec.reference_state<1||...
   spec.reference_state~=fix(spec.reference_state)||spec.reference_state>numel(ensemble.spin_config)
    error('cce_interlace:options','A valid reference_state index is required.');
end
if ~isfield(spec,'max_states')||~isscalar(spec.max_states)||~isfinite(spec.max_states)||...
   spec.max_states<1||spec.max_states~=fix(spec.max_states)
    error('cce_interlace:budget','An explicit positive integer state budget is required.');
end
reference=ensemble.spin_config(spec.reference_state);
if ~isequal(reference.bath_ids,parameters.cce.bath_spins)||...
   numel(reference.assignment)~=numel(parameters.cce.bath_spins)
    error('cce:spin_config','Interlaced reference bath IDs and assignments must match the spin_system.');
end
for j=1:numel(parameters.cce.bath_spins)
    dimension=spin_system.comp.mults(parameters.cce.bath_spins(j));
    index=reference.assignment(j);
    if ~isfinite(index)||index~=fix(index)||index<1||index>dimension
        error('cce:spin_config','Interlaced reference assignments must select local basis states.');
    end
    vector=ensemble.bases{j}(:,index);
    if ~isequal(size(reference.bath_rho{j}),[dimension dimension])||...
       norm(reference.bath_rho{j}-vector*vector','fro')>dimension^2*eps
        error('cce:spin_config',...
              'Interlaced reference bath_rho must agree with its basis assignment; use explicit_product for mixed states.');
    end
end
end


