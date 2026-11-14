% Compile central and product bath states in stable physical-ID order.
% Syntax:
%
%                    ensemble=cce_states(model,parameters)
%
% bath_state.mode is mixed_product, explicit_product, sampled, interlaced.
% Explicit product states use bath_state.rho cells in model.bath_ids order.
% Sampling supplies basis/probabilities cells (or the character 'uniform'),
% samples and seed. Threefry substreams are global state IDs, not workers.
% Each context contains central_rho, bath_ids, bath_rho, expectations,
% state_id, configuration_id and weight. Exterior fields use these same
% states for every cluster. No global random stream is modified.
%
% ilya.kuprov@weizmann.ac.il

function ensemble=cce_states(model,parameters)

% Validate explicit physical states and sampling controls
grumble(model,parameters);

% Resolve the central state without introducing a hidden polarization
d0=prod(model.multiplicities(model.central_ids));
if isfield(parameters.central,'rho')
    central=parameters.central.rho;
else
    state=(parameters.central.alpha+parameters.central.beta)/sqrt(2);
    central=state*state';
end

% Compile local densities or explicitly selected sampling distributions
ids=model.bath_ids; nb=numel(ids); mode=parameters.bath_state.mode;
rhos=cell(1,nb); bases=cell(1,nb); probabilities=cell(1,nb);
for j=1:nb
    d=model.multiplicities(ids(j));
    switch mode
        case 'mixed_product'
            rhos{j}=eye(d)/d;
        case 'explicit_product'
            rhos{j}=parameters.bath_state.rho{j};
        otherwise
            if ischar(parameters.bath_state.basis)
                bases{j}=eye(d);
            else
                bases{j}=parameters.bath_state.basis{j};
            end
            if ischar(parameters.bath_state.probabilities)
                probabilities{j}=ones(d,1)/d;
            else
                probabilities{j}=parameters.bath_state.probabilities{j}(:);
            end
    end
end

% Bind deterministic identities before any parallel work can be scheduled
ns=1; config_id=1;
if ismember(mode,{'sampled','interlaced'}), ns=parameters.bath_state.samples; end
if isfield(parameters.bath_state,'configuration_id')
    config_id=parameters.bath_state.configuration_id;
end
contexts=repmat(struct('central_rho',central,'bath_ids',ids,...
                       'bath_rho',{rhos},'expectations',zeros(nb,3),...
                       'state_id',0,'configuration_id',config_id,'weight',1/ns),1,ns);
assignments=zeros(ns,nb);
for s=1:ns

    % Assign the same global product state to every cluster in this sample
    if ismember(mode,{'sampled','interlaced'})
        stream=RandStream('Threefry','Seed',parameters.bath_state.seed);
        stream.Substream=(config_id-1)*ns+s;
        for j=1:nb
            index=find(rand(stream)<=cumsum(probabilities{j}),1);
            vector=bases{j}(:,index); rhos{j}=vector*vector';
            assignments(s,j)=index;
        end
    end
    contexts(s).bath_rho=rhos; contexts(s).state_id=s;

    % Calculate exterior spin expectations in the local Cartesian basis
    for j=1:nb
        d=model.multiplicities(ids(j)); operators=local_spin(d);
        for a=1:3
            contexts(s).expectations(j,a)=real(trace(rhos{j}*operators{a}));
        end
    end
end

% Return explicit ensemble records rather than an implicit random iterator
ensemble=struct('central_rho',central,'central_dimension',d0,...
                'bath_ids',ids,'contexts',contexts,'weights',[contexts.weight],...
                'assignments',assignments,'mode',mode,...
                'configuration_id',config_id,'state_ids',1:ns);

end

% Construct standard descending-m spin matrices without Spinach embedding
function operators=local_spin(d)
s=(d-1)/2; m=s:-1:-s;
raising=diag(sqrt(s*(s+1)-m(2:end).*(m(2:end)+1)),1);
operators={(raising+raising')/2,(raising-raising')/(2i),diag(m)};
end

% Check a finite positive trace-one local density matrix
function check_density(rho,d)
if ~isnumeric(rho)||~isequal(size(rho),[d d])||any(~isfinite(rho),'all')
    error('Density dimensions or finite entries are invalid.');
end
bound=d^2*eps*max(1,norm(rho,'fro'));
if norm(rho-rho','fro')>bound||abs(trace(rho)-1)>bound||...
   min(eig((rho+rho')/2))<-bound
    error('Density must be Hermitian, positive semidefinite and trace one.');
end
end

% Validate the full state specification before generating assignments
function grumble(model,parameters)
if ~isstruct(model)||~all(isfield(model,{'central_ids','bath_ids','multiplicities'}))||...
   ~isstruct(parameters)||~isfield(parameters,'central')||...
   ~isfield(parameters,'bath_state')||~isfield(parameters.bath_state,'mode')
    error('Canonical model, central state and bath_state.mode are required.');
end
known={'mode','rho','probabilities','basis','samples','seed','configuration_id','temperature'};
if ~isempty(setdiff(fieldnames(parameters.bath_state),known))
    error('Unknown bath_state field.');
end
if isfield(parameters.bath_state,'temperature')&&~isempty(parameters.bath_state.temperature)
    error('cce_states:capability','Thermal preparation requires explicit local densities.');
end
mode=parameters.bath_state.mode;
if ~ismember(mode,{'mixed_product','explicit_product','sampled','interlaced'})
    error('Unknown bath state mode.');
end
d=prod(model.multiplicities(model.central_ids));
if isfield(parameters.central,'rho')
    check_density(parameters.central.rho,d);
elseif ~all(isfield(parameters.central,{'alpha','beta'}))||...
       ~strcmp(parameters.method,'cce')
    error('central.rho is required except for conditional branch coherence.');
else
    cce_project(model,parameters);
end
nb=numel(model.bath_ids);
if strcmp(mode,'explicit_product')
    if ~isfield(parameters.bath_state,'rho')||~iscell(parameters.bath_state.rho)||...
       numel(parameters.bath_state.rho)~=nb
        error('bath_state.rho must contain one local density per bath site.');
    end
    for j=1:nb
        check_density(parameters.bath_state.rho{j},model.multiplicities(model.bath_ids(j)));
    end
end
if ismember(mode,{'sampled','interlaced'})
    if ~all(isfield(parameters.bath_state,{'samples','seed','basis','probabilities'}))
        error('Sampling requires samples, seed, basis and probabilities.');
    end
    for name={'samples','seed'}
        value=parameters.bath_state.(name{1});
        if ~isnumeric(value)||~isscalar(value)||~isfinite(value)||value<0||value~=fix(value)
            error('Sampling counts and seed must be nonnegative integers.');
        end
    end
    if parameters.bath_state.samples<1||parameters.bath_state.seed>2^32-1
        error('samples must be positive and seed must fit uint32.');
    end
    for name={'basis','probabilities'}
        value=parameters.bath_state.(name{1});
        if ~(ischar(value)&&strcmp(value,'uniform'))&&...
           ~(iscell(value)&&numel(value)==nb)
            error('Sampling distributions must be cells or explicit uniform.');
        end
    end
    for j=1:nb
        d=model.multiplicities(model.bath_ids(j));
        if iscell(parameters.bath_state.basis)
            b=parameters.bath_state.basis{j};
            if ~isnumeric(b)||~isequal(size(b),[d d])||any(~isfinite(b),'all')||...
               norm(b'*b-eye(d),'fro')>d^2*eps
                error('Each sampling basis must be orthonormal.');
            end
        end
        if iscell(parameters.bath_state.probabilities)
            p=parameters.bath_state.probabilities{j};
            if ~isnumeric(p)||~isreal(p)||numel(p)~=d||any(~isfinite(p))||...
               any(p<0)||abs(sum(p)-1)>d*eps
                error('Sampling probabilities must be normalized and nonnegative.');
            end
        end
    end
end
if isfield(parameters.bath_state,'configuration_id')
    id=parameters.bath_state.configuration_id;
    if ~isnumeric(id)||~isreal(id)||~isscalar(id)||~isfinite(id)||id<1||id~=fix(id)
        error('configuration_id must be a positive integer.');
    end
end
end

