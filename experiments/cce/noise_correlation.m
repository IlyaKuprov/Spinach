% Bath fluctuation autocorrelation with explicit quantum ordering. Syntax:
%
%       [A,aux]=noise_correlation(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    spin_system - native retained subsystem from cce_subsystem
%
%    parameters  - experiment structure with the fields described below
%
%    H,R,K       - angular-frequency generators in the retained basis;
%                  L=H+1i*R+1i*K for conventional matrix generators
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
% squared operator units. aux always retains moment, mean_time, and mean_zero;
% ensemble means must be formed before connected subtraction. Fourier
% transformation is separate: cce_spectrum uses exp(+1i*omega*t), omega in rad/s.
%
% Outputs:
%
%    A   - correlation column in squared operator units, with one entry
%          per requested time
%
%    aux - raw moments, first moments, ordering, validity, and absolute
%          floating-point error estimates; connect after ensemble averaging
%
% ilya.kuprov@weizmann.ac.il

function [A,aux]=noise_correlation(spin_system,parameters,H,R,K)

% Require explicit ordering and physical fluctuation operator semantics
grumble(spin_system,parameters,H,R,K);

% Get the retained dimensions and physical spin identifiers
noise=parameters.noise;
obs_dim=spin_system.observed_dimension;
bath_dim=spin_system.bath_dimension;
dimension=obs_dim*bath_dim;
ids=spin_system.physical_ids;

% Form the product bath density and the retained fluctuation operator
B=zeros(dimension);
density=1;
for j=1:numel(spin_system.bath_rho)
    density=kron(density,spin_system.bath_rho{j});
end
if isfield(noise,'terms')
    for term=reshape(noise.terms,1,[])

        % Embed each supported physical fluctuation term
        if all(ismember(term.support,ids))
            B=B+embed(term.operator,term.support,ids,spin_system.multiplicities);
        end

    end
else
    local=spin_system;
    local.formalism='zeeman-hilb';
    for descriptor=reshape(noise.descriptors,1,[])

        % Restrict the named spin sum to this retained cluster
        selected=descriptor.spins(ismember(descriptor.spins,ids));
        if ~isempty(selected)
            B=B+descriptor.coefficient*...
                cce_operator(local,descriptor.operators,selected);
        end

    end
end
B=B(1:bath_dim,1:bath_dim);

% Construct the declared bath Hamiltonian without tracing an observed transition
system=spin_system;
system.formalism='zeeman-hilb';
local=parameters;
local.cce.approximation='generalized';
if isfield(local,'bath_jumps')
    local=rmfield(local,'bath_jumps');
end
[ham,~,~]=cce_generators(system,local,[]);
if strcmp(noise.manifold,'bath_only')
    ham=reshape(full(ham),[bath_dim obs_dim bath_dim obs_dim]);
    bath_ham=zeros(bath_dim);
    for n=1:obs_dim
        bath_ham=bath_ham+reshape(ham(:,n,:,n),bath_dim,bath_dim)/obs_dim;
    end
else
    transition=system.model.transition;
    if strcmp(noise.manifold,'alpha')
        vector=transition.alpha;
    else
        vector=transition.beta;
    end
    projector=kron(vector,eye(bath_dim));
    bath_ham=projector'*ham*projector;
    if norm(ham*projector-projector*bath_ham,'fro')>dimension*eps*max(1,norm(ham,'fro'))
        error('noise_correlation:manifold','The fixed noise manifold is not invariant.');
    end
end

% Restrict bath-only jump dynamics to the selected observed-spin block
loss=zeros(bath_dim^2);
if ~isempty(R)&&any(R,'all')
    slots=reshape(1:dimension^2,dimension,dimension);
    slots=slots(1:bath_dim,1:bath_dim);
    loss=R(slots(:),slots(:));
end

% Allocate the raw moments and arithmetic error bounds
n_times=numel(parameters.times);
moment=complex(zeros(n_times,1));
mean_time=moment;
mean_zero=trace(B*density);
bound=zeros(n_times,1);
if ~any(loss,'all')

    % Preserve the Hilbert-space arithmetic and its unitary round-off estimate
    [vectors,energies]=eig((bath_ham+bath_ham')/2,'vector');
    for n=1:n_times

        % Evaluate the ordered or symmetrised Heisenberg moment
        U=vectors*diag(exp(-1i*energies*parameters.times(n)))*vectors';
        evolved=U'*B*U;
        moment(n)=trace(evolved*B*density);
        if strcmp(noise.ordering,'symmetrized')
            moment(n)=(moment(n)+trace(B*evolved*density))/2;
        end
        mean_time(n)=trace(evolved*density);

    end
    bound=64*eps*bath_dim*(1+norm(bath_ham,2)*parameters.times(:))*norm(B,2)^2;

else

    % Dissipative correlations use the full bath Liouville generator
    L=kron(eye(bath_dim),bath_ham)-kron(bath_ham.',eye(bath_dim))+1i*loss;
    initial=B*density;
    if strcmp(noise.ordering,'symmetrized')
        initial=(initial+density*B)/2;
    end

    % Scale source formation and trace contraction in the Hilbert-Schmidt norm
    roundoff=64*eps*numel(B);
    detector=norm(B,'fro');
    purity_norm=norm(density,'fro');
    source_norm=norm(initial,'fro');
    initial_error=roundoff*detector*purity_norm;
    frequency=norm(L,1);
    for n=1:n_times

        % Propagate the correlation source and the bath density
        time=parameters.times(n);
        propagator=expm(-1i*L*time);
        moment(n)=trace(B*reshape(propagator*initial(:),bath_dim,bath_dim));
        mean_time(n)=trace(B*reshape(propagator*density(:),bath_dim,bath_dim));

        % Transport source uncertainty through the actual nonunitary action
        gain=sqrt(norm(propagator,1)*norm(propagator,inf));
        arithmetic=roundoff*gain*(1+(time>0)+time*frequency);
        bound(n)=detector*(gain*initial_error+arithmetic*source_norm);
        if noise.connected

            % Both ensemble first moments are at most the detector norm in magnitude
            mean_error=arithmetic*detector*purity_norm;
            bound(n)=bound(n)+detector*(mean_error+initial_error)+...
                     mean_error*(roundoff*detector)+eps*(abs(moment(n))+2*detector^2);

        end

    end

end

% Preserve the raw moments before any connected subtraction
A=moment;
if noise.connected
    A=A-mean_time*mean_zero;
end
aux=struct('moment',moment,'mean_time',mean_time,'mean_zero',mean_zero,...
           'ordering',noise.ordering,'connected',noise.connected,...
           'valid',isfinite(A),'abs_error_bound',bound);

end

% Embed a possibly nonconsecutive ordered support in the retained tensor product
function A=embed(local,support,ids,dims)

% Tensor the local operator with the complementary identity
[~,positions]=ismember(support,ids);
rest=setdiff(1:numel(ids),positions,'stable');
order=[positions rest];
A=kron(local,eye(prod(dims(rest))));
if numel(ids)>1

    % Restore the physical spin order on both tensor indices
    tensor=reshape(A,[fliplr(dims(order)) fliplr(dims(order))]);
    inverse=zeros(size(order));
    inverse(order)=1:numel(order);
    permutation=numel(ids)+1-fliplr(inverse);
    A=reshape(permute(tensor,[permutation numel(ids)+permutation]),...
              prod(dims),prod(dims));

end

end

% Reject unspecified connected conventions and non-Hermitian observables
function grumble(system,parameters,H,R,K)
if isfield(parameters,'sequence')
    error('noise_correlation:sequence',...
          'Use finite_drive, refocusing_train, or hahn_echo for controls.');
end
if ~isstruct(system)||~isfield(system,'kind')||~strcmp(system.kind,'cce_subsystem')||...
   ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)||~isfield(parameters,'noise')||...
   ~all(isfield(parameters.noise,{'ordering','connected','manifold'}))||...
   (isfield(parameters.noise,'terms')==isfield(parameters.noise,'descriptors'))
    error('noise_correlation:input',...
          'A native subsystem and complete noise convention are required.');
end
noise=parameters.noise;
if ~ismember(noise.ordering,{'ordered','symmetrized'})||...
   ~ismember(noise.manifold,{'bath_only','alpha','beta'})||...
   ~islogical(noise.connected)||~isscalar(noise.connected)
    error('noise_correlation:convention',...
          'Declare ordering, connected flag and bath manifold.');
end
if ~isfield(parameters,'times')||any(~isfinite(parameters.times))||...
   any(parameters.times<0)||any(diff(parameters.times)<0)
    error('noise_correlation:times',...
          'Elapsed times must be finite, nonnegative and increasing.');
end
dimension=system.observed_dimension*system.bath_dimension;
if 16*(8*dimension^2+8*system.bath_dimension^4+4*numel(parameters.times))>...
   parameters.cce.numerics.max_memory_bytes
    error('noise_correlation:memory',...
          'Noise propagation exceeds the declared memory budget.');
end
if isfield(noise,'terms')
    for term=noise.terms
        dimension=prod(system.model.multiplicities(term.support));
        if ~all(ismember(term.support,system.model.bath_ids))||...
           ~isequal(size(term.operator),[dimension dimension])||...
           any(~isfinite(term.operator),'all')||...
           norm(term.operator-term.operator','fro')>100*eps*norm(term.operator,'fro')
            error('noise_correlation:operator',...
                  'Fluctuation operators must be Hermitian and supported on bath spins.');
        end
        overlap=ismember(term.support,system.physical_ids);
        if any(overlap)&&~all(overlap)
            error('noise_correlation:support',...
                  'A fluctuation support is split by this cluster.');
        end
    end
else
    descriptors=parameters.noise.descriptors;
    if ~isstruct(descriptors)||isempty(descriptors)||~isvector(descriptors)||...
       ~all(isfield(descriptors,{'operators','spins','coefficient'}))
        error('noise_correlation:descriptors',...
              'Supply a nonempty array of complete noise descriptors.');
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
            error('noise_correlation:descriptors',...
                  'Use Hermitian spin labels, unique bath IDs and a finite real coefficient.');
        end
    end
end
if isfield(parameters,'bath_jumps')
    for jump=reshape(parameters.bath_jumps,1,[])
        if ~all(ismember(jump.support,system.model.bath_ids))
            error('noise_correlation:jumps',...
                  'Noise propagation requires bath-only jumps.');
        end
    end
end
if ~strcmp(noise.manifold,'bath_only')&&~isfield(system.model,'transition')
    error('noise_correlation:transition',...
          'A selected noise manifold requires an observed transition.');
end
end


