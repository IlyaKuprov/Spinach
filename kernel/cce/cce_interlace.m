% Conditional product-state averages for interlaced CCE (W12). Syntax:
%
%       raw=cce_interlace(model,cluster,ensemble,schedule,parameters)
%
% parameters.interlaced contains max_states (explicit enumeration budget),
% reference_state (ensemble context index), and optionally condition_sites
% (cell array of physical-ID sets). Otherwise K is the union of cylinders
% associated with included superclusters of cluster. Denominator subsets
% MUST reuse raw.condition_sites, not their own default cylinder family.
% The returned raw value is the probability-weighted conditional average;
% configurations are unique and weights are conditional probabilities.
% The distribution signature records assignments and weights exactly.
%
% Inputs:
%   model, cluster, ensemble, schedule, parameters - physical model, selected
%   cluster, complete sampled ensemble, compiled sequence and interlaced
%   settings.
%   parameters.interlaced - max_states, reference_state (ensemble context
%   index), optional condition_sites cell array; an optional evaluator
%   function handle supports open-system/noise dispatch.
%
% Outputs:
%   raw - conditional raw cluster observable with condition_sites and
%   conditional-average diagnostics for consistent denominator reduction.
%
% See also cce_states, cce_reduce.
% Reference: https://doi.org/10.1103/PhysRevB.86.035452
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_interlace(model,cluster,ensemble,schedule,parameters)

% Validate the explicit conditioning and local product distribution
grumble(model,cluster,ensemble,parameters);

% Locate included superclusters unless a shared denominator K was supplied
spec=parameters.interlaced; reference=ensemble.contexts(spec.reference_state);
reference.assignment=ensemble.assignments(spec.reference_state,:);
bases=cell(size(model.bath_ids)); probabilities=bases;
for j=1:numel(model.bath_ids)
    d=model.multiplicities(model.bath_ids(j));
    if ischar(parameters.bath_state.basis), bases{j}=eye(d); else, bases{j}=parameters.bath_state.basis{j}; end
    if ischar(parameters.bath_state.probabilities)
        probabilities{j}=ones(d,1)/d;
    else
        probabilities{j}=parameters.bath_state.probabilities{j}(:);
    end
end
if isfield(spec,'condition_sites')
    cylinders=spec.condition_sites;
else
    family=cce_clusters(model,cce_graph(model,parameters),parameters);
    cylinders={family(arrayfun(@(c) all(ismember(cluster.sites,c.sites)),family)).sites};
end
if isempty(cylinders), error('cce_interlace:empty','No conditioning supercluster.'); end

% Remove contained cylinders without changing the conditioning union
keep=true(size(cylinders));
for j=1:numel(cylinders)
    for k=1:numel(cylinders)
        if j~=k&&all(ismember(cylinders{j},cylinders{k}))&&...
           (numel(cylinders{j})<numel(cylinders{k})||j>k)
            keep(j)=false; break;
        end
    end
end
cylinders=cylinders(keep);

% Reject invalid cylinder IDs and bound enumeration before allocating it
counts=zeros(size(cylinders)); dimensions=model.multiplicities(model.bath_ids);
for k=1:numel(cylinders)
    if ~isnumeric(cylinders{k})||~isreal(cylinders{k})||numel(unique(cylinders{k}))~=numel(cylinders{k})||...
       any(~ismember(cylinders{k},model.bath_ids))||~all(ismember(cluster.sites,cylinders{k}))
        error('cce_interlace:cylinder','Invalid conditioning supercluster.');
    end
    counts(k)=prod(dimensions(ismember(model.bath_ids,cylinders{k})));
end
if sum(counts)>spec.max_states
    error('cce_interlace:budget','Cylinder enumeration exceeds interlaced.max_states.');
end

% Enumerate each cylinder around the same global assignment
assignments=zeros(sum(counts),numel(model.bath_ids)); cursor=0;
for k=1:numel(cylinders)
    positions=find(ismember(model.bath_ids,cylinders{k}));
    block=repmat(reference.assignment,counts(k),1); radix=1;
    for j=positions
        block(:,j)=mod(floor((0:counts(k)-1)'/radix),dimensions(j))+1;
        radix=radix*dimensions(j);
    end
    assignments(cursor+(1:counts(k)),:)=block; cursor=cursor+counts(k);
end

% Deduplicate overlapping cylinders before computing physical probabilities
assignments=unique(assignments,'rows','sorted');
logweights=zeros(size(assignments,1),1);
for j=1:numel(model.bath_ids)
    logweights=logweights+log(probabilities{j}(assignments(:,j)));
end
shift=max(logweights);
if ~isfinite(shift), error('cce_interlace:zero_mass','Conditioning event has zero probability.'); end
weights=exp(logweights-shift); logmass=shift+log(sum(weights));
weights=weights/sum(weights); mass=exp(logmass);

% Evaluate each distinct configuration once using the unchanged cluster engine
raw=[];
for k=1:numel(weights)
    if weights(k)==0, continue; end
    context=reference; context.assignment=assignments(k,:);
    for j=1:numel(model.bath_ids)
        vector=bases{j}(:,assignments(k,j));
        context.bath_rho{j}=vector*vector';
        spin=pauli(dimensions(j)); x=spin.x; y=spin.y; z=spin.z; operators={x,y,z};
        for a=1:3
            context.expectations(j,a)=real(trace(context.bath_rho{j}*operators{a}));
        end
    end
    item=cce_cluster(model,cluster,context,schedule,parameters);
    if isempty(raw)
        raw=item; raw.value=zeros(size(item.value));
        raw.abs_error_bound=zeros(size(item.abs_error_bound));
        raw.density=zeros(size(item.density));
    end
    raw.value=raw.value+weights(k)*item.value;
    raw.abs_error_bound=raw.abs_error_bound+weights(k)*item.abs_error_bound;
    raw.density=raw.density+weights(k)*item.density;
    raw.valid=raw.valid&item.valid;
end

% Form the irreducible W12 factor using precisely the same conditioning event
if (~isfield(spec,'raw_only')||~spec.raw_only)&&~isempty(cluster.sites)
    family=cce_clusters(model,cce_graph(model,parameters),parameters);
    subsets=family(arrayfun(@(entry) all(ismember(entry.sites,cluster.sites)),family));
    factors=cell(1,numel(subsets)); errors=factors; masks=factors;
    conditioned=parameters; conditioned.interlaced.condition_sites=cylinders;
    conditioned.interlaced.raw_only=true;
    empty=cce_interlace(model,struct('id',0,'sites',[]),ensemble,schedule,conditioned);
    for j=1:numel(subsets)
        if isequal(subsets(j).sites,cluster.sites)
            item=raw;
        else
            item=cce_interlace(model,subsets(j),ensemble,schedule,conditioned);
        end
        denominator=empty.value; uncertainty=empty.abs_error_bound; valid=empty.valid;
        for k=1:j-1
            if all(ismember(subsets(k).sites,subsets(j).sites))
                uncertainty=abs(denominator).*errors{k}+abs(factors{k}).*uncertainty+uncertainty.*errors{k};
                denominator=denominator.*factors{k}; valid=valid&masks{k};
            end
        end
        resolved=valid&item.valid&abs(denominator)>uncertainty+parameters.numerics.atol+...
                 parameters.numerics.rtol*abs(denominator);
        factors{j}=NaN(size(item.value)); errors{j}=Inf(size(item.value)); masks{j}=resolved;
        factors{j}(resolved)=item.value(resolved)./denominator(resolved);
        errors{j}(resolved)=(item.abs_error_bound(resolved)+...
            abs(factors{j}(resolved)).*uncertainty(resolved))./...
            (abs(denominator(resolved))-uncertainty(resolved));
    end
    target=find(arrayfun(@(entry) isequal(sort(entry.sites),sort(cluster.sites)),subsets));
    if numel(target)~=1, error('cce_interlace:family','Requested cluster absent from included family.'); end
    raw.factor=factors{target}; raw.factor_valid=masks{target}&isfinite(raw.factor);
    raw.factor_abs_error_bound=errors{target};
end

% Return the exact conditioning signature for safe reuse by denominator subsets
raw.condition_sites=cylinders; raw.assignments=assignments;
raw.weights=weights; raw.condition_mass=mass; raw.log_condition_mass=logmass;
raw.distribution_signature=struct('assignments',assignments,'weights',weights);

end

% Check all public inputs before enumerating any states
function grumble(model,cluster,ensemble,parameters)
if ~isstruct(model)||~all(isfield(model,{'bath_ids','multiplicities'}))||...
   ~isstruct(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   ~isstruct(ensemble)||~all(isfield(ensemble,{'contexts','assignments'}))||...
   ~isstruct(parameters)||~isfield(parameters,'interlaced')||...
   ~isfield(parameters,'bath_state')||~all(isfield(parameters.bath_state,{'basis','probabilities'}))
    error('cce_interlace:input','Canonical model, cluster, ensemble and interlaced options required.');
end
spec=parameters.interlaced;
if ~isstruct(spec)||~all(isfield(spec,{'reference_state','max_states'}))||...
   ~isscalar(spec.reference_state)||~isreal(spec.reference_state)||...
   ~isfinite(spec.reference_state)||spec.reference_state<1||...
   spec.reference_state~=fix(spec.reference_state)||spec.reference_state>numel(ensemble.contexts)||...
   ~isscalar(spec.max_states)||~isreal(spec.max_states)||~isfinite(spec.max_states)||...
   spec.max_states<1||spec.max_states~=fix(spec.max_states)
    error('cce_interlace:options','Finite positive integer reference_state and max_states required.');
end
if isfield(spec,'raw_only')&&(~islogical(spec.raw_only)||~isscalar(spec.raw_only))
    error('cce_interlace:options','raw_only must be logical.');
end
if isfield(spec,'condition_sites')&&~iscell(spec.condition_sites)
    error('cce_interlace:options','condition_sites must be a cell array.');
end
assignment=ensemble.assignments(spec.reference_state,:);
if numel(assignment)~=numel(model.bath_ids)||any(assignment<1)||...
   any(assignment>model.multiplicities(model.bath_ids))||any(assignment~=fix(assignment))
    error('cce_interlace:assignment','Interlacing requires a sampled product-state assignment.');
end
end


