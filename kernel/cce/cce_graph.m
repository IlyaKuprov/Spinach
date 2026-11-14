% Physical bath connectivity and deterministic selection scores. Syntax:
%
%                  graph=cce_graph(model,parameters)
%
% model contains physical_ids, central_ids, bath_ids, multiplicities,
% coordinates_A, and pairs: records ids,tensor, where ids is a two-element physical
% ID row and tensor is a real Cartesian coupling in rad/s. Empty pairs
% must be an empty struct array. Coordinates are in Angstrom.
%
% parameters.graph.mode is complete, distance, coupling, or explicit;
% radius, threshold, adjacency specify the latter three respectively.
% score is coupling_norm (default), resonance (mismatch matrix and
% positive mismatch_scale, both rad/s), or supplied (scores matrix).
% bath_radius defaults to Inf. Hamiltonian cutoffs never enter here.
%
% graph contains bath_ids, sparse logical adjacency, scores, pair_norm,
% central_mediated, approximate, omitted_bath_ids, selection_mode, and
% requires_convergence. Matrix ordering
% follows bath_ids. gCCE retains edges between bath sites attached to
% the same connected central register unless approximate=true.
%
% Inputs:
%   model - physical IDs, central/bath IDs, multiplicities, coordinates in
%   Angstrom and pair records (ids, real Cartesian tensor in rad/s).
%   parameters.graph - mode; mode-specific radius (Angstrom), threshold
%   (rad/s), or adjacency; score is coupling_norm (default), resonance
%   (mismatch matrix and positive mismatch_scale, both rad/s), or supplied
%   (scores matrix); optional bath_radius defaults to Inf.
%
% Outputs:
%   graph - bath_ids, sparse logical adjacency, scores, pair_norm,
%   central_mediated, approximate, omitted_bath_ids, selection_mode and
%   requires_convergence.
%
% ilya.kuprov@weizmann.ac.il

function graph=cce_graph(model,parameters)

% Validate inputs
grumble(model,parameters);

% Locate physical sites without interpreting IDs as array subscripts
[~,bath_pos]=ismember(model.bath_ids,model.physical_ids);
[~,central_pos]=ismember(model.central_ids,model.physical_ids);
ids=model.bath_ids; n=numel(ids); nc=numel(central_pos);

% Obtain physical two-spin operator norms before graph pruning
pair_norm=zeros(n); attachments=false(n,nc); central_graph=eye(nc)>0;
for k=1:numel(model.pairs)

    % Resolve this unordered physical pair
    pair=model.pairs(k);
    [~,a]=ismember(pair.ids(1),model.physical_ids);
    [~,b]=ismember(pair.ids(2),model.physical_ids);
    ma=model.multiplicities(a); mb=model.multiplicities(b);
    spins=cell(2,3);
    for q=1:2

        % Construct angular momentum matrices at physical multiplicity
        d=[ma mb]; d=d(q); s=(d-1)/2; m=(s:-1:-s).';
        plus=diag(sqrt(s*(s+1)-m(2:end).*(m(2:end)+1)),1);
        spins(q,:)={(plus+plus')/2,(plus-plus')/(2i),diag(m)};

    end

    % Sum Cartesian terms and measure the Hermitian operator norm
    interaction=zeros(ma*mb);
    for u=1:3
        for v=1:3
            interaction=interaction+pair.tensor(u,v)*kron(spins{1,u},spins{2,v});
        end
    end
    strength=norm(interaction,2);

    % Accumulate bath edges and exact central attachment topology
    ia=find(ids==pair.ids(1)); ib=find(ids==pair.ids(2));
    ca=find(model.central_ids==pair.ids(1)); cb=find(model.central_ids==pair.ids(2));
    if ~isempty(ia)&&~isempty(ib)
        pair_norm(ia,ib)=strength; pair_norm(ib,ia)=strength;
    elseif ~isempty(ia)&&~isempty(cb)
        attachments(ia,cb)=strength>0;
    elseif ~isempty(ib)&&~isempty(ca)
        attachments(ib,ca)=strength>0;
    elseif ~isempty(ca)&&~isempty(cb)&&strength>0
        central_graph(ca,cb)=true; central_graph(cb,ca)=true;
    end

end

% Connect bath sites mediated by any connected central component
for k=1:nc
    central_graph=central_graph|(central_graph(:,k)&central_graph(k,:));
end
mediated=(double(attachments)*double(central_graph)*double(attachments.'))>0;
mediated(1:n+1:end)=false;

% Construct scores without modifying the Hamiltonian or edge norms
scores=pair_norm;
if isfield(parameters.graph,'score')
    switch parameters.graph.score
        case 'resonance'
            scores=pair_norm./hypot(1,parameters.graph.mismatch/...
                                     parameters.graph.mismatch_scale);
        case 'supplied'
            scores=parameters.graph.scores;
    end
end
scores(1:n+1:end)=0;

% Apply only the explicitly requested graph criterion
switch parameters.graph.mode
    case 'complete'
        adjacency=~eye(n);
    case 'distance'
        xyz=model.coordinates_A(bath_pos,:);
        adjacency=false(n);
        for k=1:n
            adjacency(k,:)=vecnorm(xyz-xyz(k,:),2,2)<=parameters.graph.radius;
        end
    case 'coupling'
        adjacency=(scores>0)&(scores>=parameters.graph.threshold);
    case 'explicit'
        adjacency=full(parameters.graph.adjacency);
end
adjacency(1:n+1:end)=false;

% Preserve mediated correlations unless graph approximation is explicit
approximate=isfield(parameters.graph,'approximate')&&parameters.graph.approximate;
if ismember(parameters.method,{'gcce','megcce'})&&~approximate
    adjacency=adjacency|mediated;
end

% Restrict candidate bath by distance from the union of central sites
keep=true(1,n);
if isfield(parameters,'bath_radius')&&isfinite(parameters.bath_radius)
    distance=inf(n,1);
    for k=central_pos
        distance=min(distance,vecnorm(model.coordinates_A(bath_pos,:)-...
                                     model.coordinates_A(k,:),2,2));
    end
    keep=distance.'<=parameters.bath_radius;
end

% Return physical ordering and explicit approximation diagnostics
graph=struct('bath_ids',ids(keep),'adjacency',sparse(adjacency(keep,keep)),...
             'scores',scores(keep,keep),'pair_norm',pair_norm(keep,keep),...
             'central_mediated',sparse(mediated(keep,keep)),...
             'approximate',approximate,'omitted_bath_ids',ids(~keep),...
             'selection_mode',parameters.graph.mode,...
             'requires_convergence',~strcmp(parameters.graph.mode,'complete')||any(~keep));

end

% Validate shapes, physical IDs, graph semantics, and physical units
function grumble(model,parameters)
if ~isstruct(model)||~isscalar(model)||~all(isfield(model,{'physical_ids',...
        'central_ids','bath_ids','multiplicities','coordinates_A','pairs'}))
    error('cce_graph:model','incomplete scalar model structure.');
end
ids=model.physical_ids; n=numel(model.bath_ids);
if ~isnumeric(ids)||~isreal(ids)||(~isrow(ids)&&~isempty(ids))||...
   any(~isfinite(ids))||any(ids<1)||any(ids~=fix(ids))||numel(unique(ids))~=numel(ids)
    error('cce_graph:ids','physical IDs must be distinct positive integers.');
end
for name={'central_ids','bath_ids'}
    values=model.(name{1});
    if ~isnumeric(values)||~isreal(values)||(~isrow(values)&&~isempty(values))||...
       numel(unique(values))~=numel(values)||~all(ismember(values,ids))
        error('cce_graph:ids','central and bath IDs must be unique physical ID rows.');
    end
end
if ~isempty(intersect(model.central_ids,model.bath_ids))||isempty(model.central_ids)
    error('cce_graph:ids','central IDs must be nonempty and disjoint from bath IDs.');
end
if ~isnumeric(model.multiplicities)||~isreal(model.multiplicities)||...
   numel(model.multiplicities)~=numel(ids)||any(~isfinite(model.multiplicities))||...
   any(model.multiplicities<1)||any(model.multiplicities~=fix(model.multiplicities))
    error('cce_graph:dimensions','multiplicities must be positive integers.');
end
if ~isnumeric(model.coordinates_A)||~isreal(model.coordinates_A)||...
   ~isequal(size(model.coordinates_A),[numel(ids) 3])||any(isinf(model.coordinates_A),'all')
    error('cce_graph:coordinates','coordinates_A must be a physical-site by three array.');
end
if ~isstruct(model.pairs)||(~isempty(model.pairs)&&~all(isfield(model.pairs,{'ids','tensor'})))
    error('cce_graph:pairs','pairs must be records ids,tensor.');
end
pairs=zeros(numel(model.pairs),2);
for k=1:numel(model.pairs)
    p=model.pairs(k);
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
if ~isstruct(parameters)||~isscalar(parameters)||~all(isfield(parameters,{'method','graph'}))||...
   ~ismember(parameters.method,{'cce','gcce','mecce','megcce','noise'})||...
   ~isstruct(parameters.graph)||~isscalar(parameters.graph)||...
   ~isfield(parameters.graph,'mode')||...
   ~ismember(parameters.graph.mode,{'complete','distance','coupling','explicit'})
    error('cce_graph:parameters','explicit method and graph mode are required.');
end
g=parameters.graph;
for field={'radius','threshold'}
    if (strcmp(g.mode,'distance')&&strcmp(field{1},'radius'))||...
       (strcmp(g.mode,'coupling')&&strcmp(field{1},'threshold'))
        if ~isfield(g,field{1})||~isnumeric(g.(field{1}))||~isreal(g.(field{1}))||...
           ~isscalar(g.(field{1}))||isnan(g.(field{1}))||g.(field{1})<0
            error('cce_graph:cutoff','nonnegative graph radius or threshold required.');
        end
    end
end
if strcmp(g.mode,'explicit')&&(~isfield(g,'adjacency')||~islogical(g.adjacency)||...
   ~issparse(g.adjacency)||~isequal(size(g.adjacency),[n n])||...
   ~isequal(g.adjacency,g.adjacency.')||any(diag(g.adjacency)))
    error('cce_graph:adjacency','adjacency must be symmetric sparse logical with zero diagonal.');
end
if isfield(g,'approximate')&&(~islogical(g.approximate)||~isscalar(g.approximate))
    error('cce_graph:approximate','approximate must be scalar logical.');
end
if isfield(g,'score')
    if ~ismember(g.score,{'coupling_norm','resonance','supplied'})
        error('cce_graph:score','unknown graph score.');
    end
    if strcmp(g.score,'resonance')
        if ~all(isfield(g,{'mismatch','mismatch_scale'}))||...
           ~isnumeric(g.mismatch_scale)||~isreal(g.mismatch_scale)||...
           ~isscalar(g.mismatch_scale)||~isfinite(g.mismatch_scale)||g.mismatch_scale<=0||...
           ~isnumeric(g.mismatch)||~isreal(g.mismatch)||~isequal(size(g.mismatch),[n n])||...
           any(~isfinite(g.mismatch),'all')||~isequal(abs(g.mismatch),abs(g.mismatch.'))
            error('cce_graph:score','resonance needs symmetric absolute mismatches and positive scale.');
        end
    elseif strcmp(g.score,'supplied')
        if ~isfield(g,'scores')||~isnumeric(g.scores)||~isreal(g.scores)||...
           ~isequal(size(g.scores),[n n])||any(~isfinite(g.scores),'all')||...
           any(g.scores<0,'all')||~isequal(g.scores,g.scores.')
            error('cce_graph:score','supplied scores must be finite symmetric nonnegative.');
        end
    end
end
if isfield(parameters,'bath_radius')&&(~isnumeric(parameters.bath_radius)||...
   ~isreal(parameters.bath_radius)||~isscalar(parameters.bath_radius)||...
   isnan(parameters.bath_radius)||parameters.bath_radius<0)
    error('cce_graph:radius','bath_radius must be nonnegative.');
end
if (strcmp(g.mode,'distance')||(isfield(parameters,'bath_radius')&&...
    isfinite(parameters.bath_radius)))&&any(~isfinite(model.coordinates_A),'all')
    error('cce_graph:coordinates','distance selection requires finite physical coordinates.');
end
end
