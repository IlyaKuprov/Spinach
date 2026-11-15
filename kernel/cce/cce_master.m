% Exact Lindblad propagation of one ME-CCE or ME-gCCE cluster. Syntax:
%
%          raw=cce_master(model,cluster,spin_config,schedule,parameters)
%
% Interfaces match cce_cluster; jumps follow cce_liouvillian. Supports FID,
% ideal physical pulses; mecce controls must preserve or swap manifolds.
% Finite controls require the controls extension and are refused here.
% Static exterior coherent fields are retained, not exterior dissipators.
% Generalized initial states and normalization use spin_config.rho0.
% Multiplicative populations are not a generally valid T1 approximation.
%
% Inputs:
%   model, cluster, spin_config, schedule, parameters - local physical/state
%   inputs, method, ideal sequence, explicit dissipation jumps and numerical
%   limits.
%
% Outputs:
%   raw - cluster observable/density, validity and numerical error metadata
%   for cce_reduce.
%
% See also cce_liouvillian, cce_observe.
% Reference: https://doi.org/10.1103/PhysRevLett.132.250401
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_master(model,cluster,spin_config,schedule,parameters)

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters,model,cluster.sites,spin_config);

% Validate the propagation contract before allocating cluster matrices
grumble(model,cluster,spin_config,schedule,parameters);
sites=cluster.sites; hc=struct('backend','auto');

% Retain the same static exterior assignment as closed-system clusters
if strcmp(parameters.meanfield,'static')
    mask=~ismember(spin_config.bath_ids,sites); exterior=spin_config.bath_ids(mask);
    if isfield(parameters,'meanfield_radius')&&isfinite(parameters.meanfield_radius)
        distance=inf(1,numel(exterior));
        for id=model.observed_ids
            delta=model.coordinates_A(exterior,:)-model.coordinates_A(id,:);
            distance=min(distance,sqrt(sum(delta.^2,2)).');
        end
        exterior=exterior(distance<=parameters.meanfield_radius);
    end
    [~,where]=ismember(exterior,spin_config.bath_ids);
    hc.exterior_ids=exterior; hc.exterior_expectations=spin_config.expectations(where,:);
end

% Construct the native generator and the cluster product bath state
[G,meta]=cce_liouvillian(model,sites,hc,parameters);
db=meta.bath_dimension; d_obs=meta.observed_dimension; rho_b=1;
for site=sites
    where=find(spin_config.bath_ids==site,1); rho_b=kron(rho_b,spin_config.bath_rho{where});
end
nt=numel(schedule.times); densities=complex(zeros(d_obs,d_obs,nt)); values=complex(zeros(nt,1));

% Select the vector exponential-action manifold without rebuilding a basis
engine=model.spin_system; engine.bas.formalism='sphten-liouv';

% Disable stock entry chopping because the CCE bound assumes full precision
engine.tols.prop_chop=0;

% Compile ideal physical controls before selecting the propagated block
pulses=cell(1,numel(schedule.events));
for j=1:numel(pulses), pulses{j}=pulse_matrix(schedule.events(j),meta.hamiltonian,model); end
if meta.projected
    rho_initial=rho_b;
    flip=reshape(1:db^2,db,db).'; flip=flip(:);
    alternate=conj(G(flip,flip)); va=kron(parameters.transition.alpha,eye(db));
    vb=kron(parameters.transition.beta,eye(db)); blocks=cell(1,numel(pulses));
    for j=1:numel(pulses)
        P=pulses{j}; aa=va'*P*va; ab=va'*P*vb; ba=vb'*P*va; bb=vb'*P*vb;
        tolerance=100*eps*meta.dimension;
        if norm(ab,'fro')+norm(ba,'fro')<tolerance
            blocks{j}={aa,bb,false};
        elseif norm(aa,'fro')+norm(bb,'fro')<tolerance
            blocks{j}={ba,ab,true};
        else
            error('cce:manifold_mixing','The pulse mixes the selected spin manifolds; use megcce.');
        end
        pair=blocks{j};
        if norm(pair{1}'*pair{1}-eye(db),'fro')+...
           norm(pair{2}'*pair{2}-eye(db),'fro')>tolerance
            error('cce:manifold_mixing','The pulse mixes the selected spin manifolds; use megcce.');
        end
    end
else
    rho_initial=kron(spin_config.rho0,rho_b);
end

% Evolve chronologically with the existing Spinach propagation utilities
for k=1:nt
    vector=rho_initial(:); previous=0; manifold=false;
    if schedule.times(k)>0
        for j=schedule.order{k}(:).'
            dt=schedule.positions(j,k)-previous;
            if meta.projected&&manifold, generator=alternate; else, generator=G; end
            vector=evolve(engine,generator,vector,dt,parameters);
            if meta.projected
                pair=blocks{j};
                if manifold, pair(1:2)=pair([2 1]); end
                density=reshape(vector,db,db); density=pair{1}*density*pair{2}';
                vector=density(:); manifold=xor(manifold,pair{3});
            else
                density=reshape(vector,meta.dimension,meta.dimension);
                density=pulses{j}*density*pulses{j}'; vector=density(:);
            end
            previous=schedule.positions(j,k);
        end
    end
    if meta.projected&&manifold, generator=alternate; else, generator=G; end
    vector=evolve(engine,generator,vector,schedule.times(k)-previous,parameters);
    if meta.projected
        values(k)=trace(reshape(vector,db,db));
    else
        density=reshape(vector,meta.dimension,meta.dimension);
        for a=1:d_obs
            for b=1:d_obs
                densities(a,b,k)=trace(density((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
    end
end

% Do not expose an uncomputed observed density for projected propagation
if meta.projected, densities=[]; end

% Observe using the supplied initial state as the normalization reference
parameters.rho0=spin_config.rho0;
if ~meta.projected
    values=cce_observe(densities,parameters);
end

% Bound arithmetic error explicitly without claiming a physical CCE bound
scale=ones(1,size(values,2));
if strcmp(parameters.quantity,'observable')
    operators=parameters.coil;
    if ~iscell(operators), operators={operators}; end
    scale=cellfun(@(op) norm(op,2),operators(:).');
end
if ~meta.projected&&isfield(parameters,'normalized')&&parameters.normalized
    observe=parameters; observe.normalized=false;
    initial=cce_observe(parameters.rho0,observe);
    scale=scale./abs(initial);
end
argument=norm(G,1)*schedule.times(:);
bound=(64*eps*size(G,1)*(1+argument))*scale;
raw=struct('cluster_id',cluster.id,'value',values,...
           'valid',isfinite(values),'abs_error_bound',bound,...
           'density',densities,'generator',meta,...
           'error_bound_kind','floating_point_estimate');

end

% Enforce supported controls and a valid spin configuration
function grumble(model,cluster,spin_config,schedule,parameters)
if ~isstruct(model)||~isstruct(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   ~isstruct(spin_config)||~all(isfield(spin_config,{'rho0','bath_ids','bath_rho','expectations'}))||...
   ~isstruct(schedule)||~all(isfield(schedule,{'times','events','positions'}))
    error('Prepared model, cluster, bath spin_config and schedule are required.');
end
if ~isstruct(parameters)||~all(isfield(parameters,{'method','quantity','meanfield'}))||...
   ~ismember(parameters.method,{'mecce','megcce'})
    error('Master propagation requires an explicit mecce or megcce method.');
end
if strcmp(parameters.method,'mecce')&&~strcmp(parameters.quantity,'signal')
    error('ME-CCE returns the selected-transition signal only.');
end
if isfield(parameters.numerics,'propagator')&&...
   ~strcmp(parameters.numerics.propagator,'stock')
    error('ME propagation accepts stock; unspecified uses step exp-action.');
end
if strcmp(parameters.meanfield,'static')
    for event=schedule.events
        if any(ismember(event.target,setdiff(spin_config.bath_ids,cluster.sites)))
            error('cce:capability','Exterior bath pulses require time-dependent mean fields.');
        end
    end
end
end

% Apply a native-generator exponential with an explicit stock dense option
function vector=evolve(engine,G,vector,time,parameters)
if isfield(parameters.numerics,'propagator')&&strcmp(parameters.numerics.propagator,'stock')
    vector=propagator(engine,1i*G,time)*vector;
else
    vector=step(engine,1i*G,vector,time);
end
end

% Form a physical ideal pulse in the existing observed-first tensor order
function pulse=pulse_matrix(event,meta,model)
d=meta.dimension; generator=zeros(d);
if isfield(event,'operator')&&~isempty(event.operator)

    % Scale the generator by the optional rotation angle
    if isfield(event,'angle')&&~isempty(event.angle), event.operator=event.angle*event.operator; end
    if isequal(event.target,model.observed_ids)
        generator=kron(event.operator,eye(meta.bath_dimension));
    elseif isscalar(event.target)
        position=find(meta.physical_ids==event.target,1);
        if ~isempty(position)
            left=prod(meta.dimensions(1:position-1));
            right=prod(meta.dimensions(position+1:end));
            generator=kron(eye(left),kron(event.operator,eye(right)));
        end
    else
        error('cce_master:operator_target','Operator pulses must target the whole observed-spin system or one physical site.');
    end
else
    axis=event.axis;
    if ischar(axis), axis=double('xyz'==axis); end
    for id=event.target
        position=find(meta.physical_ids==id,1);
        if isempty(position), continue; end
        for a=1:3
            generator=generator+event.angle*axis(a)*meta.operators{position,a};
        end
    end
end
pulse=expm(-1i*generator);
end


