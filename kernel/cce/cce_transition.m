% Resolve an ordered observed-spin transition in the actual laboratory frame.
%
% Syntax:
%
%     [transition,meta]=cce_transition(spin_system,parameters,assumptions)
%
% Parameters:
%
%   spin_system - canonical create output or a compiled cce_model; tensors
%                 already use rad/s
%   parameters  - observed is a row of physical spin IDs; transition.spins
%                 is the same ordered row. Choose exactly one selector:
%                 ms is n-observed-by-2, levels is a 1-by-2 row of distinct
%                 ascending-energy indices, or alpha/beta are explicit
%                 orthonormal columns in descending-m product ordering.
%                 orientation is a 1-by-3 Euler-angle row in radians
%                 (the zero-angle frame is used when it is absent).
%   assumptions - stock assume name or cell row of assume arguments
%
% Outputs:
%
%   transition  - ordered alpha/beta columns; spins preserves physical IDs
%   meta        - selector, ordering, energies, signed gap, and nearest gaps
%                 in rad/s, eigenstate residuals, and the selected frame
%
% Only the observed subsystem is diagonalised. Energy-index selection is
% rejected when either state belongs to an unresolved degenerate manifold;
% explicit orthonormal alpha/beta states then specify the physical choice.
% Non-eigenstate selectors have NaN nearest gaps and explicit residuals.
% No optical population or nuclear polarisation is prepared here.
%
% ilya.kuprov@weizmann.ac.il

function [transition,meta]=cce_transition(spin_system,parameters,assumptions)

% Check the selection before constructing the observed subsystem
grumble(spin_system,parameters,assumptions);

% A compiled spin system reuses the retained generator; ordinary callers use stock
model=[];
if isfield(spin_system,'cce')
    model=spin_system;
end

% Preserve the caller's random stream across stock subsystem construction
rng_state=rng;
rng_restore=onCleanup(@()rng(rng_state));
selection=parameters.transition;
ids=parameters.observed;
dims=spin_system.comp.mults(ids);
dimension=prod(dims);

% Construct only the observed Hamiltonian in the declared physical frame
angles=[0 0 0];
if isfield(parameters,'orientation')
    angles=parameters.orientation;
end
if ischar(assumptions)
    assumption_args={assumptions};
else
    assumption_args=assumptions;
end
native=~isempty(model)&&isfield(parameters.transition,'ms')&&...
       strcmp(model.bas.formalism,'zeeman-hilb')&&...
       ~any(isfield(parameters,{'bath_state','bath_samples'}))&&...
       (~isfield(parameters,'bath_jumps')||isempty(parameters.bath_jumps));
if ~native
    system=spin_system;
    if isfield(system,'bas')
        system=rmfield(system,'bas');
    end
    system=kill_spin(system,setdiff(1:spin_system.comp.nspins,ids));
    bas.formalism='zeeman-hilb';
    bas.approximation={'none'};
    system=basis(system,bas);
    system=assume(system,assumption_args{:});
    if isfield(parameters,'orientation')
        angles=parameters.orientation;
    else
        angles=[0 0 0];
    end
    [H,Q]=hamiltonian(system);
    H=full(H+orientation(Q,angles));

    % Preserve the caller's tensor-product ordering rather than sorting spin IDs
    if numel(ids)>1
        [~,order]=ismember(ids,sort(ids));
        indices=reshape(1:dimension,fliplr(spin_system.comp.mults(sort(ids))));
        indices=permute(indices,numel(ids)+1-fliplr(order));
        indices=indices(:);
        H=H(indices,indices);
    end

else

    % Reuse the compiled observed subsystem without exterior mean fields
    p=parameters;
    p.cce.meanfield='none';
    p.cce.approximation='generalized';
    model.bas.formalism='zeeman-hilb';
    config=struct('configuration_id',1,'state_id',1,'weight',1);
    [system,~]=cce_subsystem(model,p,[],config,[]);
    [H,~,~,~]=cce_generators(system,p,[]);
    H=full(H);
end

% Check Hermiticity at the rounding scale of the observed Hamiltonian
tolerance=dimension*eps(max(1,norm(H,'fro')));
if norm(H-H','fro')>tolerance
    error('cce_transition:hermiticity','The observed Hamiltonian must be Hermitian.');
end
[V,D]=eig((H+H')/2,'vector');
[energies,order]=sort(real(D));
V=V(:,order);

% Resolve product magnetic quantum numbers or energy-ordered eigenstates
if isfield(selection,'ms')
    vectors=ones(1,2);
    for n=1:numel(ids)
        columns=zeros(dims(n),2);
        for k=1:2
            index=(dims(n)-1)/2-selection.ms(n,k)+1;
            columns(index,k)=1;
        end
        vectors=[kron(vectors(:,1),columns(:,1)) ...
                 kron(vectors(:,2),columns(:,2))];
    end
    kind='ms';
elseif isfield(selection,'levels')
    for k=selection.levels
        others=setdiff(1:dimension,k);
        if any(abs(energies(k)-energies(others))<=tolerance)
            error('cce_transition:degeneracy',...
                  'Energy indices do not resolve this manifold; supply explicit alpha/beta states.');
        end
    end
    vectors=V(:,selection.levels);
    for k=1:2
        [~,pivot]=max(abs(vectors(:,k)));
        vectors(:,k)=vectors(:,k)*conj(vectors(pivot,k))/abs(vectors(pivot,k));
    end
    kind='levels';
else
    vectors=[selection.alpha selection.beta];
    kind='explicit';
end

% Return the ordered states and their physically meaningful energy metadata
transition=struct('spins',ids,'alpha',vectors(:,1),'beta',vectors(:,2));
meta=struct('selector',kind,'spins',ids,'multiplicities',dims,...
            'ordering','descending_m','orientation',angles,...
            'assumptions',{assumption_args},'energy_units','rad/s');
meta.energies=real(diag(vectors'*H*vectors)).';
meta.gap=meta.energies(2)-meta.energies(1);
meta.residuals=zeros(1,2);
meta.nearest_gaps=NaN(1,2);
meta.eigenvalues=energies.';
for k=1:2
    meta.residuals(k)=norm(H*vectors(:,k)-meta.energies(k)*vectors(:,k));
    if meta.residuals(k)<=tolerance
        gaps=sort(abs(energies-meta.energies(k)));
        meta.nearest_gaps(k)=gaps(2);
    end
end
meta.degeneracy_tolerance=tolerance;

end

% Validate physical spin IDs and exactly one descriptive transition selector
function grumble(spin_system,p,assumptions)
if ~isstruct(p)||~isscalar(p)||~all(isfield(p,{'observed','transition'}))
    error('cce_transition:inputs','parameters.observed and transition are required.');
end
ids=p.observed;
if ~isnumeric(ids)||~isreal(ids)||~isrow(ids)||isempty(ids)||...
   any(~isfinite(ids))||any(ids~=fix(ids))||any(ids<1)||...
   any(ids>spin_system.comp.nspins)||numel(unique(ids))~=numel(ids)
    error('cce_transition:spins','Observed spin IDs must be a distinct valid row.');
end
s=p.transition;
if ~isstruct(s)||~isscalar(s)||~isfield(s,'spins')||~isequal(s.spins,ids)||...
   ~isempty(setdiff(fieldnames(s),{'spins','ms','levels','alpha','beta'}))
    error('cce_transition:selector','The selector must name the ordered observed spins and supported fields.');
end
explicit=isfield(s,'alpha')||isfield(s,'beta');
if double(isfield(s,'ms'))+double(isfield(s,'levels'))+double(explicit)~=1
    error('cce_transition:selector','Choose exactly one of ms, levels or alpha/beta.');
end
dims=spin_system.comp.mults(ids);
dimension=prod(dims);
if isfield(s,'ms')
    if ~isnumeric(s.ms)||~isreal(s.ms)||~isequal(size(s.ms),[numel(ids) 2])||...
       any(~isfinite(s.ms),'all')
        error('cce_transition:ms','ms must be a finite n-observed-by-2 real array.');
    end
    indices=(dims(:)-1)/2-s.ms+1;
    if any(indices<1,'all')||any(indices>dims(:),'all')||...
       any(indices~=fix(indices),'all')||isequal(s.ms(:,1),s.ms(:,2))
        error('cce_transition:ms','Invalid or identical magnetic-quantum-number states.');
    end
elseif isfield(s,'levels')
    if ~isnumeric(s.levels)||~isreal(s.levels)||~isequal(size(s.levels),[1 2])||...
       any(~isfinite(s.levels))||any(s.levels~=fix(s.levels))||...
       any(s.levels<1)||any(s.levels>dimension)||s.levels(1)==s.levels(2)
        error('cce_transition:levels','levels must be two distinct valid energy indices.');
    end
else
    if ~all(isfield(s,{'alpha','beta'}))||~isnumeric(s.alpha)||~isnumeric(s.beta)||...
       ~isequal(size(s.alpha),[dimension 1])||~isequal(size(s.beta),[dimension 1])||...
       any(~isfinite(s.alpha))||any(~isfinite(s.beta))
        error('cce_transition:states','Explicit states must be finite columns of the observed dimension.');
    end
    vectors=[s.alpha s.beta];
    if norm(vectors'*vectors-eye(2),'fro')>dimension*eps*max(1,norm(vectors,'fro')^2)
        error('cce_transition:states','Explicit transition states must be orthonormal.');
    end
end
if isfield(p,'orientation')&&(~isnumeric(p.orientation)||~isreal(p.orientation)||...
   ~isequal(size(p.orientation),[1 3])||any(~isfinite(p.orientation)))
    error('cce_transition:orientation','orientation must be a finite Euler-angle row in radians.');
end
if ~(ischar(assumptions)&&isrow(assumptions)&&~isempty(assumptions))&&...
   ~(iscell(assumptions)&&isrow(assumptions)&&~isempty(assumptions)&&ischar(assumptions{1}))
    error('cce_transition:assumptions','Supply the explicit stock assume arguments.');
end
end


