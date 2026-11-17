% Bath fluctuation autocorrelation with explicit quantum ordering. Syntax:
%
%       [A,aux]=noise_correlation(spin_system,parameters,H,R,K)
%
% parameters.noise contains Hermitian terms (support in physical spin IDs,
% operator in that support order), ordering ('ordered' or 'symmetrized'),
% connected (logical), and manifold ('bath_only', 'alpha', or 'beta').
% Instead of terms, noise.descriptors is a nonempty structure array with
% operators ('Lx', 'Ly', 'Lz', or 'E'), spins (a nonempty physical bath-ID
% row), and coefficient (finite real scalar, in the fluctuation units).
% Each descriptor is coefficient times the sum of the selected one-spin
% operators. Descriptors on exterior spins contribute zero in that cluster;
% the full bath-ID selection is validated before restriction. Do not mix
% descriptors with explicit terms. No full-bath operator is constructed.
% parameters.times contains elapsed times in seconds. A is a column in
% squared operator units. aux always retains moment, mean_time and mean_zero;
% ensemble means must be formed before connected subtraction. Fourier
% transformation is separate: cce_spectrum uses exp(+1i*omega*t), omega in rad/s.
%
% ilya.kuprov@weizmann.ac.il

function [A,aux]=noise_correlation(spin_system,parameters,H,R,K)

% Require explicit ordering and physical fluctuation operator semantics
grumble(spin_system,parameters,H,R,K);
noise=parameters.noise; dc=spin_system.observed_dimension;
db=spin_system.bath_dimension; d=dc*db; ids=spin_system.physical_ids;
bytes=16*(8*d^2+8*db^4+4*numel(parameters.times));
if bytes>parameters.cce.numerics.max_memory_bytes
    error('noise_correlation:memory','Noise propagation exceeds the declared memory budget.');
end
B=zeros(d); density=1;
for j=1:numel(spin_system.bath_rho), density=kron(density,spin_system.bath_rho{j}); end
if isfield(noise,'terms')
    for term=reshape(noise.terms,1,[])
        overlap=ismember(term.support,ids);
        if ~any(overlap), continue; end
        if ~all(overlap), error('noise_correlation:support','A fluctuation support is split by this cluster.'); end
        B=B+embed(term.operator,term.support,ids,spin_system.multiplicities);
    end
else
    local=spin_system; local.formalism='zeeman-hilb';
    for descriptor=reshape(noise.descriptors,1,[])
        selected=descriptor.spins(ismember(descriptor.spins,ids));
        if ~isempty(selected)
            B=B+descriptor.coefficient*cce_operator(local,descriptor.operators,selected);
        end
    end
end
B=B(1:db,1:db);

% Construct the declared bath Hamiltonian without tracing an observed transition
system=spin_system; system.formalism='zeeman-hilb';
local=parameters; local.cce.approximation='generalized';
if isfield(local,'bath_jumps'), local=rmfield(local,'bath_jumps'); end
[ham,~,~]=cce_generators(system,local,[]);
if strcmp(noise.manifold,'bath_only')
    ham=reshape(full(ham),[db dc db dc]); hb=zeros(db);
    for n=1:dc, hb=hb+reshape(ham(:,n,:,n),db,db)/dc; end
else
    if ~isfield(system.model,'transition')
        error('noise_correlation:transition','A selected noise manifold requires an observed transition.');
    end
    transition=system.model.transition;
    if strcmp(noise.manifold,'alpha'), vector=transition.alpha; else, vector=transition.beta; end
    projector=kron(vector,eye(db)); hb=projector'*ham*projector;
    if norm(ham*projector-projector*hb,'fro')>d*eps*max(1,norm(ham,'fro'))
        error('noise_correlation:manifold','The fixed noise manifold is not invariant.');
    end
end

% Restrict bath-only jump dynamics to the selected observed-spin block
loss=zeros(db^2);
if ~isempty(R)&&any(R,'all')
    slots=reshape(1:d^2,d,d); slots=slots(1:db,1:db); loss=R(slots(:),slots(:));
end
nt=numel(parameters.times); moment=complex(zeros(nt,1)); mean_time=moment;
mean_zero=trace(B*density); bound=zeros(nt,1);
if ~any(loss,'all')

    % Preserve the Hilbert-space arithmetic and its unitary round-off estimate
    [vectors,energies]=eig((hb+hb')/2,'vector');
    for n=1:nt
        U=vectors*diag(exp(-1i*energies*parameters.times(n)))*vectors';
        evolved=U'*B*U; moment(n)=trace(evolved*B*density);
        if strcmp(noise.ordering,'symmetrized')
            moment(n)=(moment(n)+trace(B*evolved*density))/2;
        end
        mean_time(n)=trace(evolved*density);
    end
    bound=64*eps*db*(1+norm(hb,2)*parameters.times(:))*norm(B,2)^2;

else

    % Dissipative correlations use the full bath Liouville generator
    L=kron(eye(db),hb)-kron(hb.',eye(db))+1i*loss;
    initial=B*density;
    if strcmp(noise.ordering,'symmetrized'), initial=(initial+density*B)/2; end
    for n=1:nt
        propagator=expm(-1i*L*parameters.times(n));
        moment(n)=trace(B*reshape(propagator*initial(:),db,db));
        mean_time(n)=trace(B*reshape(propagator*density(:),db,db));
    end

end
A=moment;
if noise.connected, A=A-mean_time*mean_zero; end
aux=struct('moment',moment,'mean_time',mean_time,'mean_zero',mean_zero,...
           'ordering',noise.ordering,'connected',noise.connected,...
           'valid',isfinite(A),'abs_error_bound',bound);

end

% Embed a possibly nonconsecutive ordered support in the retained tensor product
function A=embed(local,support,ids,dims)
[~,positions]=ismember(support,ids); rest=setdiff(1:numel(ids),positions,'stable');
order=[positions rest]; A=kron(local,eye(prod(dims(rest))));
if numel(ids)>1
    tensor=reshape(A,[fliplr(dims(order)) fliplr(dims(order))]);
    inverse=zeros(size(order)); inverse(order)=1:numel(order);
    permutation=numel(ids)+1-fliplr(inverse);
    A=reshape(permute(tensor,[permutation numel(ids)+permutation]),prod(dims),prod(dims));
end
end

% Reject unspecified connected conventions and non-Hermitian observables
function grumble(system,p,H,R,K)
if ~isstruct(system)||~isfield(system,'kind')||~strcmp(system.kind,'cce_subsystem')||...
   ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)||~isfield(p,'noise')||...
   ~all(isfield(p.noise,{'ordering','connected','manifold'}))||...
   (isfield(p.noise,'terms')==isfield(p.noise,'descriptors'))
    error('noise_correlation:input','A native subsystem and complete noise convention are required.');
end
noise=p.noise;
if ~ismember(noise.ordering,{'ordered','symmetrized'})||...
   ~ismember(noise.manifold,{'bath_only','alpha','beta'})||...
   ~islogical(noise.connected)||~isscalar(noise.connected)
    error('noise_correlation:convention','Declare ordering, connected flag and bath manifold.');
end
if ~isfield(p,'times')||any(~isfinite(p.times))||any(p.times<0)||any(diff(p.times)<0)
    error('noise_correlation:times','Elapsed times must be finite, nonnegative and increasing.');
end
if isfield(noise,'terms')
    for term=noise.terms
        d=prod(system.model.multiplicities(term.support));
        if ~all(ismember(term.support,system.model.bath_ids))||...
           ~isequal(size(term.operator),[d d])||any(~isfinite(term.operator),'all')||...
           norm(term.operator-term.operator','fro')>100*eps*norm(term.operator,'fro')
            error('noise_correlation:operator','Fluctuation operators must be Hermitian and supported on bath spins.');
        end
    end
else
    descriptors=p.noise.descriptors;
    if ~isstruct(descriptors)||isempty(descriptors)||~isvector(descriptors)||...
       ~all(isfield(descriptors,{'operators','spins','coefficient'}))
        error('noise_correlation:descriptors','Supply a nonempty array of complete noise descriptors.');
    end
    for descriptor=reshape(descriptors,1,[])
        if ~ischar(descriptor.operators)||~isrow(descriptor.operators)||...
           ~ismember(descriptor.operators,{'Lx','Ly','Lz','E'})||...
           ~isnumeric(descriptor.spins)||~isreal(descriptor.spins)||...
           ~isrow(descriptor.spins)||isempty(descriptor.spins)||...
           any(~ismember(descriptor.spins,system.model.bath_ids))||...
           numel(unique(descriptor.spins))~=numel(descriptor.spins)||...
           ~isnumeric(descriptor.coefficient)||~isreal(descriptor.coefficient)||...
           ~isscalar(descriptor.coefficient)||~isfinite(descriptor.coefficient)
            error('noise_correlation:descriptors','Use Hermitian spin labels, unique bath IDs and a finite real coefficient.');
        end
    end
end
if isfield(p,'bath_jumps')
    for jump=reshape(p.bath_jumps,1,[])
        if ~all(ismember(jump.support,system.model.bath_ids))
            error('noise_correlation:jumps','Noise propagation requires bath-only jumps.');
        end
    end
end
end


