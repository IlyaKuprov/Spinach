% Projected product-state averages for interlaced CCE (W12). Syntax:
%
%       raw=cce_interlace(model,cluster,ensemble,schedule,parameters,evaluator)
%
% parameters.interlaced contains max_states (explicit enumeration budget),
% reference_state (ensemble spin_config index), and optionally condition_sites
% (cell array of physical-ID sets). Otherwise K is the union of cylinders
% associated with included superclusters of cluster. Denominator subsets
% MUST reuse raw.condition_sites, not their own default cylinder family.
% The reference bath densities must agree with their basis assignment.
% The returned raw value is the probability-weighted conditional average;
% configurations are unique and weights are conditional probabilities.
% The distribution signature records assignments and weights exactly.
%
% Inputs:
%   model, cluster, ensemble, schedule, parameters - physical model, selected
%   cluster, complete sampled ensemble, compiled sequence and interlaced
%   settings.
%   parameters.interlaced - max_states, reference_state (ensemble spin_config
%   index), optional condition_sites cell array and raw_only flag. The
%   declared method selects closed-system or master-equation propagation.
%
% Outputs:
%   raw - conditional raw cluster observable with condition_sites and
%   conditional-average diagnostics for consistent denominator reduction.
%
% See also cce_states, cce_reduce.
% Reference: https://doi.org/10.1103/PhysRevB.86.035452
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_interlace(model,cluster,ensemble,schedule,parameters,evaluator)

% Preserve an explicitly supplied experiment evaluator through every recursion
if nargin<6
    evaluator=@cce_cluster;
    if ismember(parameters.method,{'mecce','megcce'}), evaluator=@cce_master; end
end
if ~isa(evaluator,'function_handle')
    error('cce_interlace:evaluator','evaluator must be a function handle.');
end

% Native contexts carry nested algorithm controls and an explicit evaluator
if isfield(parameters,'cce')
    raw=native_interlace(model,cluster,ensemble,schedule,parameters,evaluator); return
end

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate the explicit conditioning and local product distribution
grumble(model,cluster,ensemble,parameters);

% Locate included superclusters unless a shared denominator K was supplied
spec=parameters.interlaced; reference=ensemble.spin_configs(spec.reference_state);
if isfield(reference,'assignment')&&...
   ~isequal(reference.assignment,ensemble.assignments(spec.reference_state,:))
    error('cce:spin_config','Interlaced spin_config and ensemble reference assignments must agree.');
end
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

% Validate the supplied reference before replacing any enumerated bath state
cce_validate(parameters,model,cluster.sites,reference);
if ~isequal(reference.bath_ids,model.bath_ids)||...
   numel(reference.assignment)~=numel(model.bath_ids)
    error('cce:spin_config','Interlaced reference bath IDs and assignments must match the model.');
end

% Require each reference density to represent its declared local assignment
for j=1:numel(model.bath_ids)
    dimension=model.multiplicities(model.bath_ids(j)); basis=bases{j};
    index=reference.assignment(j);
    if ~isfinite(index)||index~=fix(index)||index<1||index>dimension
        error('cce:spin_config','Interlaced reference assignments must select local basis states.');
    end
    vector=basis(:,index);
    if ~isequal(size(reference.bath_rho{j}),[dimension dimension])||...
       norm(reference.bath_rho{j}-vector*vector','fro')>dimension^2*eps
        error('cce:spin_config',...
              'Interlaced reference bath_rho must agree with its basis assignment; use explicit_product for mixed states.');
    end
end

% Resolve the cylinder family only after the reference has been validated
if isfield(spec,'condition_sites')
    cylinders=spec.condition_sites;
else
    family=cce_clusters(model,cce_graph(model,parameters),parameters);

    % Retain the reference-assignment baseline when no clusters were kept
    if isempty(family)&&isempty(cluster.sites)
        cylinders={[]};
    else
        cylinders={family(arrayfun(@(c) all(ismember(cluster.sites,c.sites)),family)).sites};
    end
end
if isempty(cylinders), error('cce_interlace:empty','No conditioning supercluster.'); end

% Bound cylinder masks, enumeration counts and index vectors before allocation
bytes=3*8*numel(cylinders);
if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
    error('cce_interlace:memory',...
          'Conditioning cylinder metadata exceed numerics.max_memory_bytes.');
end

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

% Bound assignments, a cylinder block, sorting copies and index/weight vectors
rows=sum(counts); columns=numel(model.bath_ids);
bytes=8*rows*(4*columns+5);
if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
    error('cce_interlace:memory',...
          'Conditioning assignments and sorting workspace exceed numerics.max_memory_bytes.');
end

% Enumerate each cylinder around the same global assignment
assignments=zeros(rows,columns); cursor=0;
for k=1:numel(cylinders)
    positions=find(ismember(model.bath_ids,cylinders{k}));
    block=repmat(reference.assignment,counts(k),1); radix=1;
    for j=positions
        block(:,j)=mod(floor((0:counts(k)-1)'/radix),dimensions(j))+1;
        radix=radix*dimensions(j);
    end
    assignments(cursor+(1:counts(k)),:)=block; cursor=cursor+counts(k);
end

% Release the last cylinder before sorting the conditioning union
clear block

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
    spin_config=reference; spin_config.assignment=assignments(k,:);
    for j=1:numel(model.bath_ids)
        vector=bases{j}(:,assignments(k,j));
        spin_config.bath_rho{j}=vector*vector';
        spin=pauli(dimensions(j)); x=spin.x; y=spin.y; z=spin.z; operators={x,y,z};
        for a=1:3
            spin_config.expectations(j,a)=real(trace(spin_config.bath_rho{j}*operators{a}));
        end
    end

    % Preserve the physical evaluator in every recursive conditional average
    item=evaluator(model,cluster,spin_config,schedule,parameters);
    if isempty(raw)
        raw=item; raw.value=zeros(size(item.value));
        raw.abs_error_bound=zeros(size(item.abs_error_bound));
        raw.density=zeros(size(item.density));
    end
    raw.value=raw.value+weights(k)*item.value;
    raw.abs_error_bound=raw.abs_error_bound+weights(k)*item.abs_error_bound;
    raw.density=raw.density+weights(k)*item.density;
    raw.valid=raw.valid&item.valid;

    % Preserve worst constituent residuals rather than only the first assignment
    if isfield(item,'eigen_residual')
        raw.eigen_residual=max(raw.eigen_residual,item.eigen_residual);
    end
    if isfield(item,'orthogonality_residual')
        raw.orthogonality_residual=max(raw.orthogonality_residual,item.orthogonality_residual);
    end
    if isfield(item,'unitarity_residual')
        raw.unitarity_residual=max(raw.unitarity_residual,item.unitarity_residual);
    end
end

% Form the irreducible W12 factor using precisely the same conditioning event
if (~isfield(spec,'raw_only')||~spec.raw_only)&&~isempty(cluster.sites)
    family=cce_clusters(model,cce_graph(model,parameters),parameters);
    subsets=family(arrayfun(@(entry) all(ismember(entry.sites,cluster.sites)),family));
    bytes=8*rows*(columns+1)+3*25*numel(raw.value)*numel(subsets);
    if ~isfinite(bytes)||bytes>parameters.numerics.max_memory_bytes
        error('cce_interlace:memory',...
              'Conditional subset factors exceed numerics.max_memory_bytes.');
    end
    factors=cell(1,numel(subsets)); errors=factors; masks=factors;
    conditioned=parameters; conditioned.interlaced.condition_sites=cylinders;
    conditioned.interlaced.raw_only=true;
    empty=cce_interlace(model,struct('id',0,'sites',[]),ensemble,schedule,conditioned,evaluator);
    for j=1:numel(subsets)
        if isequal(subsets(j).sites,cluster.sites)
            item=raw;
        else
            item=cce_interlace(model,subsets(j),ensemble,schedule,conditioned,evaluator);
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

% Average a raw cluster or recursively divide parent-conditioned subset records
function raw=native_interlace(model,cluster,ensemble,schedule,parameters,evaluator)
spec=parameters.cce.interlaced;
if ~isfield(spec,'max_states')||~isscalar(spec.max_states)||~isfinite(spec.max_states)||...
   spec.max_states<1||spec.max_states~=fix(spec.max_states)
    error('cce_interlace:budget','An explicit positive integer state budget is required.');
end
reference=ensemble.spin_config(spec.reference_state);

% Reject a density inconsistent with the sampled reference assignment
if ~isequal(reference.bath_ids,model.bath_ids)||...
   numel(reference.assignment)~=numel(model.bath_ids)
    error('cce:spin_config','Interlaced reference bath IDs and assignments must match the model.');
end
for j=1:numel(model.bath_ids)
    dimension=model.multiplicities(model.bath_ids(j)); index=reference.assignment(j);
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
if isfield(spec,'condition_sites')
    cylinders=spec.condition_sites;
else
    cylinders={};
    for k=1:numel(model.clusters)
        parent=model.clusters(k);
        if all(ismember(cluster.sites,parent.sites)), cylinders{end+1}=parent.sites; end %#ok<AGROW>
    end
    if isempty(cylinders), cylinders={cluster.sites}; end
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
nb=numel(model.bath_ids); dimensions=model.multiplicities(model.bath_ids);
counts=cellfun(@(sites)prod(dimensions(ismember(model.bath_ids,sites))),cylinders);
count=sum(counts);
if count>spec.max_states||8*count*max(1,nb)>parameters.cce.numerics.max_memory_bytes
    error('cce_interlace:budget','Native conditioning enumeration exceeds its explicit budget.');
end
assignments=zeros(count,nb); position=0;
for k=1:numel(cylinders)
    varied=find(ismember(model.bath_ids,cylinders{k}));
    for n=1:counts(k)
        position=position+1; assignments(position,:)=reference.assignment; index=n-1;
        for j=fliplr(varied)
            assignments(position,j)=mod(index,dimensions(j))+1; index=floor(index/dimensions(j));
        end
    end
end
assignments=unique(assignments,'rows'); count=size(assignments,1); weights=zeros(count,1);
for n=1:count
    weight=1;
    for j=1:nb, weight=weight*ensemble.probabilities{j}(assignments(n,j)); end
    weights(n)=weight;
end
if sum(weights)==0, error('cce_interlace:weight','Conditioning event has zero probability.'); end
weights=weights/sum(weights); retained=find(weights>0);
for n=reshape(retained,1,[])
    configuration=reference; configuration.assignment=assignments(n,:);
    for j=1:nb
        vector=ensemble.bases{j}(:,assignments(n,j)); configuration.bath_rho{j}=vector*vector';
    end
    item=evaluator(model,cluster,configuration,schedule,parameters);
    if n==retained(1)
        raw=item; raw.value(:)=0; raw.abs_error_bound(:)=0; raw.valid(:)=true;
    end
    raw.value=raw.value+weights(n)*item.value;
    raw.abs_error_bound=raw.abs_error_bound+weights(n)*item.abs_error_bound;
    raw.valid=raw.valid&item.valid;
end
raw.condition_sites=cylinders; raw.assignment_count=numel(retained);
raw.reference=raw.value; raw.reference_error=raw.abs_error_bound;
raw.full_bath=isequal(sort(cluster.sites),sort(model.bath_ids));
if isfield(spec,'raw_only')&&spec.raw_only, return; end

% Every denominator receives exactly the parent's conditioning event
if ~isempty(cluster.sites)&&~raw.full_bath
    local=parameters; local.cce.interlaced.condition_sites=cylinders;
    local.cce.interlaced.raw_only=true;
    family=model.clusters;
    included=arrayfun(@(member)all(ismember(member.sites,cluster.sites)),family);
    family=family(included); records=cell(1,numel(family)+1);
    empty=struct('id',0,'sites',[]);
    records{1}=native_interlace(model,empty,ensemble,schedule,local,evaluator);
    for k=1:numel(family)
        if isequal(family(k).sites,cluster.sites), records{k+1}=raw;
        else, records{k+1}=native_interlace(model,family(k),ensemble,schedule,local,evaluator); end
    end
    factors=records{1}.value; errors=records{1}.abs_error_bound;
    for k=1:numel(family)
        item=records{k+1}; denominator=records{1}.value;
        denominator_error=records{1}.abs_error_bound; valid=records{1}.valid&item.valid;
        for j=family(k).proper_subsets
            position=find([family.id]==model.clusters(j).id,1);
            if isempty(position), continue; end
            factor=records{position+1};
            denominator_error=abs(denominator).*factor.abs_error_bound+...
                              abs(factor.value).*denominator_error+...
                              denominator_error.*factor.abs_error_bound;
            denominator=denominator.*factor.value; valid=valid&factor.valid;
        end
        safe=abs(denominator)>denominator_error;
        item.abs_error_bound=(item.abs_error_bound+abs(item.value./denominator).*...
                              denominator_error)./(abs(denominator)-denominator_error);
        item.value=item.value./denominator; item.valid=valid&safe&isfinite(item.value);
        item.value(~item.valid)=NaN; item.abs_error_bound(~item.valid)=Inf;
        records{k+1}=item;
    end
    position=find([family.id]==cluster.id,1); raw=records{position+1};
    raw.reference=factors; raw.reference_error=errors;
end
raw.kind='irreducible';
if raw.full_bath, raw.kind='raw_full_bath'; end
end

function grumble(model,cluster,ensemble,parameters)
if ~isstruct(model)||~all(isfield(model,{'bath_ids','multiplicities'}))||...
   ~isstruct(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   ~isstruct(ensemble)||~all(isfield(ensemble,{'spin_configs','assignments'}))||...
   ~isstruct(parameters)||~isfield(parameters,'interlaced')||...
   ~isfield(parameters,'bath_state')||~all(isfield(parameters.bath_state,{'basis','probabilities'}))
    error('cce_interlace:input','Canonical model, cluster, ensemble and interlaced options required.');
end
spec=parameters.interlaced;
if ~isstruct(spec)||~all(isfield(spec,{'reference_state','max_states'}))||...
   ~isscalar(spec.reference_state)||~isreal(spec.reference_state)||...
   ~isfinite(spec.reference_state)||spec.reference_state<1||...
   spec.reference_state~=fix(spec.reference_state)||spec.reference_state>numel(ensemble.spin_configs)||...
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


