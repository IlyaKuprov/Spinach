% Propagate one conditional or generalized CCE cluster exactly. Syntax:
%
%          raw=cce_cluster(model,cluster,context,schedule,parameters)
%
% cluster has id and sites (empty with id=0 for the central baseline).
% context is one cce_states context. The same exterior assignment is used
% in all clusters. raw supplies the cce_reduce contract: cluster_id, value,
% valid and abs_error_bound, plus central density and spectral diagnostics.
% Conditional values are normalized branch traces Tr(Ua*rhoB*Ub').
% Generalized values come from cce_observe; full density is also retained.
% Static eigensystem propagation reports residual-derived forward bounds.
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_cluster(model,cluster,context,schedule,parameters)

% Route driven generalized clusters through the local segment compiler
if isstruct(schedule)&&isfield(schedule,'controls')&&schedule.controls
    if ~strcmp(parameters.method,'gcce')
        error('cce_cluster:controls','Driven or exterior-pulse sequences require gCCE.');
    end
    segments=cce_controls(model,sort(cluster.sites),context,schedule,parameters);
    raw=segments.raw; raw.cluster_id=uint64(cluster.id); return;
end

% Validate the fixed interfaces and supported dynamics family
grumble(model,cluster,context,schedule,parameters);

% Select a stable bath ordering and the common exterior mean field
sites=sort(cluster.sites); hcontext=struct();
if strcmp(parameters.meanfield,'static')
    mask=~ismember(context.bath_ids,sites); exterior=context.bath_ids(mask);
    if isfield(parameters,'meanfield_radius')&&isfinite(parameters.meanfield_radius)
        distance=inf(1,numel(exterior));
        for id=model.central_ids
            delta=model.coordinates_A(exterior,:)-model.coordinates_A(id,:);
            distance=min(distance,sqrt(sum(delta.^2,2)).');
        end
        keep=distance<=parameters.meanfield_radius;
        exterior=exterior(keep); indices=find(mask); mask(:)=false; mask(indices(keep))=true;
    end
    hcontext.exterior_ids=exterior;
    hcontext.exterior_expectations=context.expectations(mask,:);
end
[H,meta]=cce_hamiltonian(model,sites,hcontext);
d0=meta.central_dimension; db=meta.bath_dimension; d=meta.dimension;

% Check simultaneous dense storage before assembling states and pulses
copies=12+3*numel(meta.physical_ids)+2*numel(schedule.events);
if 16*d^2*copies>parameters.numerics.max_memory_bytes
    error('cce_cluster:memory','Simultaneous propagator storage exceeds the declared budget.');
end
rho_b=1;
for id=sites
    rho_b=kron(rho_b,context.bath_rho{find(context.bath_ids==id,1)});
end

% Compile physical pulses once, preserving the order of simultaneous events
pulses=cell(1,numel(schedule.events));
for j=1:numel(pulses)
    pulses{j}=pulse_matrix(schedule.events(j),meta,model);
end
nt=numel(schedule.times); density=complex(zeros(d0,d0,nt));
bounds=zeros(nt,1);
if strcmp(parameters.method,'cce')

    % Project the complete cluster Hamiltonian onto explicit branch states
    branches=cce_project(model,parameters);
    pa=kron(branches.alpha.state,eye(db)); pb=kron(branches.beta.state,eye(db));
    ha=pa'*H*pa; hb=pb'*H*pb;
    if branches.second_order_enabled
        ha=ha+virtual_correction(branches.alpha,model,sites,hcontext,db);
        hb=hb+virtual_correction(branches.beta,model,sites,hcontext,db);
    end
    [va,ea,ra,oa]=spectrum(ha); [vb,eb,rb,ob]=spectrum(hb);
    residual=max(ra,rb); orthogonality=max(oa,ob);
    values=complex(zeros(nt,1));

    % Require central pulses to preserve the chosen two-branch subspace
    branch_pulses=cell(size(pulses)); swaps=false(size(pulses));
    for j=1:numel(pulses)
        p=pulses{j}; aa=pa'*p*pa; ab=pa'*p*pb; ba=pb'*p*pa; bb=pb'*p*pb;
        tolerance=d^3*eps;
        if norm(aa,'fro')+norm(bb,'fro')<=tolerance
            swaps(j)=true; branch_pulses{j}={ba,ab};
        elseif norm(ab,'fro')+norm(ba,'fro')<=tolerance
            branch_pulses{j}={aa,bb};
        else
            error('cce_cluster:conditional_pulse','Conditional CCE supports only branch-preserving or branch-swapping pulses.');
        end
        pair=branch_pulses{j};
        if max(norm(pair{1}'*pair{1}-eye(db),'fro'),...
               norm(pair{2}'*pair{2}-eye(db),'fro'))>tolerance
            error('cce_cluster:conditional_leakage','Pulse leaks outside the selected central states.');
        end
    end
    for k=1:nt

        % At zero duration return the initial condition, without zero-time pulses
        ua=eye(db); ub=eye(db); branch=false; previous=0;
        if schedule.times(k)>0
            for j=reshape(schedule.order{k},1,[])
                dt=schedule.positions(j,k)-previous;
                if branch
                    ua=propagate(vb,eb,dt)*ua; ub=propagate(va,ea,dt)*ub;
                else
                    ua=propagate(va,ea,dt)*ua; ub=propagate(vb,eb,dt)*ub;
                end
                pair=branch_pulses{j};
                if branch, pair=pair([2 1]); end
                ua=pair{1}*ua; ub=pair{2}*ub;
                branch=xor(branch,swaps(j)); previous=schedule.positions(j,k);
            end
            dt=schedule.times(k)-previous;
            if branch
                ua=propagate(vb,eb,dt)*ua; ub=propagate(va,ea,dt)*ub;
            else
                ua=propagate(va,ea,dt)*ua; ub=propagate(vb,eb,dt)*ub;
            end
        end
        values(k)=trace(ua*rho_b*ub');
        bounds(k)=forward_bound(db,schedule.times(k),residual,orthogonality,...
                                meta.roundoff_bound,norm(H,'fro'),numel(pulses));
    end
else

    % Evolve the full central register and trace only the bath factors
    rho_initial=kron(context.central_rho,rho_b);
    [vectors,energies,residual,orthogonality]=spectrum(H);
    for k=1:nt
        u=eye(d); previous=0;
        if schedule.times(k)>0
            for j=reshape(schedule.order{k},1,[])
                u=pulses{j}*propagate(vectors,energies,schedule.positions(j,k)-previous)*u;
                previous=schedule.positions(j,k);
            end
            u=propagate(vectors,energies,schedule.times(k)-previous)*u;
        end
        rho=u*rho_initial*u';
        for a=1:d0
            for b=1:d0
                density(a,b,k)=trace(rho((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
        bounds(k)=forward_bound(d,schedule.times(k),residual,orthogonality,...
                                meta.roundoff_bound,norm(H,'fro'),numel(pulses));
    end
    observe_parameters=parameters; observe_parameters.central.rho=context.central_rho;
    values=cce_observe(density,observe_parameters);
    if strcmp(parameters.quantity,'observable')
        operators=parameters.operators;
        if ~iscell(operators), operators={operators}; end
        scale=cellfun(@(op) norm(op,2),operators);
    else
        scale=ones(1,size(values,2));
    end
    if isfield(parameters,'normalized')&&parameters.normalized
        observe_parameters.normalized=false;
        scale=scale./abs(cce_observe(context.central_rho,observe_parameters));
    end
    bounds=bounds.*scale;
end

% Return finite raw data with auditable propagation and Hamiltonian bounds
bounds=repmat(bounds,1,size(values,2)/size(bounds,2));
raw=struct('cluster_id',uint64(cluster.id),'value',values,...
           'valid',isfinite(values),'abs_error_bound',bounds,...
           'density',density,'sites',sites,'state_id',context.state_id,...
           'configuration_id',context.configuration_id,...
           'eigen_residual',residual,'orthogonality_residual',orthogonality,...
           'hamiltonian_bound',meta.roundoff_bound,'backend',meta.backend);

end

% Evaluate the optional central-gap second-order bath Hamiltonian
function correction=virtual_correction(branch,model,sites,context,db)
nc=numel(model.central_ids); fields=cell(1,3*nc);
for c=1:nc
    id=model.central_ids(c);
    for a=1:3
        field=zeros(db);
        for k=1:numel(sites)
            d=model.multiplicities(sites(k)); ops=spin_matrices(d);
            left=prod(model.multiplicities(sites(1:k-1)));
            right=prod(model.multiplicities(sites(k+1:end)));
            for b=1:3
                field=field+model.retained{id,sites(k)}(a,b)*...
                      kron(eye(left),kron(ops{b},eye(right)));
            end
        end
        if isfield(context,'exterior_ids')
            for k=1:numel(context.exterior_ids)
                field=field+eye(db)*(model.retained{id,context.exterior_ids(k)}(a,:)*...
                                    context.exterior_expectations(k,:).');
            end
        end
        fields{3*(c-1)+a}=field;
    end
end
correction=zeros(db);
for a=1:numel(fields)
    for b=1:numel(fields)
        correction=correction+branch.susceptibility(a,b)*fields{a}*fields{b};
    end
end
end

% Form a physical event rotation in cluster tensor-factor order
function pulse=pulse_matrix(event,meta,model)
d=meta.dimension; generator=zeros(d);
if isfield(event,'operator')&&~isempty(event.operator)
    if isequal(event.target,model.central_ids)
        generator=kron(event.operator,eye(meta.bath_dimension));
    elseif isscalar(event.target)
        position=find(meta.physical_ids==event.target,1);
        if ~isempty(position)
            left=prod(meta.dimensions(1:position-1));
            right=prod(meta.dimensions(position+1:end));
            generator=kron(eye(left),kron(event.operator,eye(right)));
        end
    else
        error('cce_cluster:operator_target','Operator events target the whole central register or one physical site.');
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

% Obtain a Hermitian eigensystem and measure its backward residuals
function [vectors,energies,residual,orthogonality]=spectrum(H)
[vectors,energies]=eig(full((H+H')/2),'vector');
residual=norm(H*vectors-vectors.*energies.','fro');
orthogonality=norm(vectors'*vectors-eye(size(H)),'fro');
end

% Apply the exact static spectral exponential
function u=propagate(vectors,energies,time)
u=(vectors.*exp(-1i*time*energies).')*vectors';
end

% Bound density/branch-trace error using residuals and arithmetic counts
function bound=forward_bound(d,time,residual,orthogonality,hbound,hnorm,npulses)
operations=(12+8*npulses)*d^3;
gamma=operations*eps/(1-operations*eps);
bound=2*time*(residual+hbound)+4*(npulses+1)*orthogonality+...
      gamma*(1+time*hnorm);
end

% Construct standard descending-m local spin operators
function operators=spin_matrices(d)
s=(d-1)/2; m=s:-1:-s;
p=diag(sqrt(s*(s+1)-m(2:end).*(m(2:end)+1)),1);
operators={(p+p')/2,(p-p')/(2i),diag(m)};
end

% Reject unsupported routes explicitly instead of silently dropping physics
function grumble(model,cluster,context,schedule,parameters)
if ~isstruct(model)||~isstruct(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   ~isstruct(context)||~all(isfield(context,{'central_rho','bath_ids','bath_rho',...
                                           'expectations','state_id','configuration_id'}))||...
   ~isstruct(schedule)||~all(isfield(schedule,{'times','events','positions','order'}))
    error('Canonical model, cluster, state context and schedule are required.');
end
if ~isstruct(parameters)||~isfield(parameters,'method')||...
   ~ismember(parameters.method,{'cce','gcce'})
    error('cce_cluster:capability','This propagator implements cce and gcce only.');
end
if ~isfield(parameters,'meanfield')||~ismember(parameters.meanfield,{'none','static'})
    error('meanfield must be explicitly none or static.');
end
if isfield(parameters.numerics,'propagator')&&~strcmp(parameters.numerics.propagator,'eig')
    error('cce_cluster:capability','Only the residual-checked eig propagator is implemented.');
end
if strcmp(parameters.method,'gcce')&&isfield(parameters,'second_order')&&...
   parameters.second_order.enabled
    error('Generalized CCE cannot add a conditional second-order correction.');
end
if strcmp(parameters.meanfield,'static')
    for event=schedule.events
        if any(ismember(event.target,setdiff(context.bath_ids,cluster.sites)))
            error('cce_cluster:exterior_pulse','Exterior bath pulses require time-dependent mean fields.');
        end
    end
end
end

