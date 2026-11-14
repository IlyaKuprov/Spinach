% Additive connected bath-noise autocorrelation for one cluster. Syntax:
%
%           raw=cce_noise(model,cluster,context,schedule,parameters)
%
% parameters.noise specifies ordering ('ordered' or 'symmetrized'),
% connected (logical), and terms: a struct array with support (ordered
% physical bath IDs) and operator (matrix in that tensor order). The
% observable is the sum of included terms; a split multi-site term is
% refused. Units are the squared units of these operators. Dynamics uses
% only the bath Hamiltonian (normalized trace over the central register),
% or its conditional projection when noise.branch is 'alpha' or 'beta'.
% Each raw is Tr[B(t)B(0)rho]-Tr[B(t)rho]Tr[B(0)rho] when connected.
% No stationarity, detailed balance, Gaussianity or T1 assumption is made.
%
% Inputs:
%   model, cluster, context, schedule, parameters - cluster inputs;
%   noise.ordering=ordered|symmetrized, logical connected, and terms with
%   ordered bath support and operator matrices; optional noise.branch selects
%   alpha or beta.
%
% Outputs:
%   raw - time-dependent autocorrelation with the standard cluster IDs,
%   validity and error fields; additive reduction is selected by method
%   noise.
%
% See also cce_reduce, cce_spectrum.
%
% ilya.kuprov@weizmann.ac.il

function raw=cce_noise(model,cluster,context,schedule,parameters)

% Validate explicit quantum ordering and the additive observable contract
grumble(model,cluster,context,schedule,parameters);
[H,~]=cce_hamiltonian(model,cluster.sites,struct('backend','auto'));
d0=prod(model.multiplicities(model.central_ids)); db=size(H,1)/d0;

% Remove the central register without retaining mediated bath terms
if isfield(parameters.noise,'branch')
    branch=cce_project(model,parameters);
    if branch.second_order_enabled
        error('cce:capability','Noise virtual-state corrections are not implemented.');
    end
    vector=kron(branch.(parameters.noise.branch).state,eye(db)); hb=vector'*H*vector;
else
    hb=zeros(db);
    for k=1:d0, hb=hb+H((k-1)*db+(1:db),(k-1)*db+(1:db))/d0; end
end

% Embed the declared bath observable in stable physical tensor order
ids=cluster.sites; dims=model.multiplicities(ids); B=zeros(db);
for k=1:numel(parameters.noise.terms)
    term=parameters.noise.terms(k); included=ismember(term.support,ids);
    if ~any(included), continue; end
    if ~all(included), error('cce:noiseSupport','Group every multi-site observable support.'); end
    positions=arrayfun(@(id) find(ids==id,1),term.support);
    rest=setdiff(1:numel(ids),positions,'stable'); perm=[positions rest];
    labels=reshape(1:db,fliplr([1 dims]));
    labels=permute(labels,numel(dims)+2-fliplr([1 perm+1])); labels=labels(:);
    local=kron(term.operator,eye(prod(dims(rest))));
    B(labels,labels)=B(labels,labels)+local;
end

% Form the same cluster product state as the dynamics expansion
rho=1;
for id=ids
    where=find(context.bath_ids==id,1); rho=kron(rho,context.bath_rho{where});
end
mean_zero=trace(B*rho); nt=numel(schedule.times); values=complex(zeros(nt,1));
mean_time=complex(zeros(nt,1)); moments=complex(zeros(nt,1));
[V,D]=eig((hb+hb')/2,'vector');

% Evaluate the explicitly ordered or symmetrized quantum correlation
for k=1:nt
    U=V*diag(exp(-1i*D*schedule.times(k)))*V'; bt=U'*B*U;
    values(k)=trace(bt*B*rho);
    if strcmp(parameters.noise.ordering,'symmetrized')
        values(k)=(values(k)+trace(B*bt*rho))/2;
    end
    mean_time(k)=trace(bt*rho); moments(k)=values(k);
    if parameters.noise.connected, values(k)=values(k)-mean_time(k)*mean_zero; end
end

% Return an additive cumulant input with a scale-sensitive roundoff estimate
bound=64*eps*db*(1+norm(hb,2)*schedule.times(:))*norm(B,2)^2;
raw=struct('cluster_id',cluster.id,'value',values,'valid',isfinite(values),...
           'abs_error_bound',bound,'ordering',parameters.noise.ordering,...
           'connected',parameters.noise.connected,'reducer','additive',...
           'noise_moment',moments,'noise_mean',mean_time,'noise_mean_zero',mean_zero);

end

% Reject ambiguous quantum correlations and unsupported controls
function grumble(model,cluster,context,schedule,parameters)
if ~isstruct(model)||~isstruct(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   ~isstruct(context)||~all(isfield(context,{'bath_ids','bath_rho'}))||...
   ~isstruct(schedule)||~all(isfield(schedule,{'times','events'}))||...
   ~isstruct(parameters)||~isfield(parameters,'noise')
    error('Prepared model, cluster, context, schedule and noise parameters are required.');
end
if ~isfield(parameters,'method')||~strcmp(parameters.method,'noise')||...
   ~isfield(parameters,'quantity')||~strcmp(parameters.quantity,'autocorrelation')
    error('cce_noise:capability','Noise evaluation requires method=noise and quantity=autocorrelation.');
end
noise=parameters.noise;
if ~isstruct(noise)||~all(isfield(noise,{'ordering','connected','terms'}))||...
   ~ismember(noise.ordering,{'ordered','symmetrized'})||...
   ~islogical(noise.connected)||~isscalar(noise.connected)||~isstruct(noise.terms)
    error('Explicit ordering, connected flag and observable terms are required.');
end
if isfield(noise,'branch')&&~ismember(noise.branch,{'alpha','beta'})
    error('noise.branch must be alpha or beta.');
end
if ~isempty(schedule.events)||~strcmp(parameters.meanfield,'none')
    error('cce:capability','Noise correlation currently requires FID and no exterior field.');
end
if isfield(parameters,'jumps')&&~isempty(parameters.jumps)
    error('cce:capability','Dissipative regression is not a closed-system noise correlation.');
end
for k=1:numel(noise.terms)
    term=noise.terms(k);
    if ~all(isfield(term,{'support','operator'}))||isempty(term.support)||...
       ~isnumeric(term.support)||~isrow(term.support)||...
       numel(unique(term.support))~=numel(term.support)||...
       any(~ismember(term.support,model.bath_ids))
        error('Noise supports must contain distinct physical bath IDs.');
    end
    n=prod(model.multiplicities(term.support));
    if ~isnumeric(term.operator)||~isequal(size(term.operator),[n n])||...
       any(~isfinite(term.operator),'all')||...
       norm(term.operator-term.operator','fro')>100*eps*norm(term.operator,'fro')
        error('Noise operators must be finite Hermitian matrices of the support dimension.');
    end
end
end


