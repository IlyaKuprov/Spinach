% Native Lindblad generator for column-vectorised density matrices. Syntax:
%
%           [G,meta]=cce_liouvillian(model,sites,context,parameters)
%
% model is returned by cce_prepare; sites lists physical bath IDs. context
% accepts cce_hamiltonian fields. parameters.method is mecce or megcce;
% parameters.jumps is a struct array with support (ordered physical IDs),
% operator (dimensionless matrix in that tensor order), and rate (s^-1).
% Evolution is exp(G*t), NOT exp(-1i*G*t). Conditional blocks require
% central-diagonal jumps; a jump split across the cluster is refused.
% Disjoint exterior jumps are omitted only without a frozen field from them;
% evolving exterior fields are refused, not replaced by a dissipator.
% meta supplies ordering, dimensions, Hamiltonians and jump records.
%
% Inputs:
%   model, sites, context, parameters - prepared cluster inputs and
%   dissipation.jumps records with ordered support IDs, operator matrix and
%   rate in s^-1.
%
% Outputs:
%   G - dense GKSL generator acting on column-vectorized density/coherence
%   blocks.
%   meta - tensor order, dimensions, Hamiltonians and included jump records.
%
% See also cce_master, cce_hamiltonian.
% Reference: https://doi.org/10.1103/PhysRevLett.132.250401
%
% ilya.kuprov@weizmann.ac.il

function [G,meta]=cce_liouvillian(model,sites,context,parameters)

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate dimensions, units and the explicit method
grumble(model,sites,context,parameters);

% Build the Hamiltonian in central-first tensor order
[H,hmeta]=cce_hamiltonian(model,sites,context);
ids=hmeta.physical_ids; dims=model.multiplicities(ids); d=size(H,1);
d0=prod(model.multiplicities(model.central_ids)); db=d/d0;
conditional=strcmp(parameters.method,'mecce');

% Guard the dense generator allocation before forming Kronecker products
if conditional, n=db^2; else, n=d^2; end
if 8*16*n^2>parameters.numerics.max_memory_bytes
    error('cce:memory','Liouvillian construction exceeds max_memory_bytes.');
end

% Assemble the coherent full-density or conditional-block generator
if conditional
    branch=cce_project(model,parameters);
    if branch.second_order_enabled
        error('cce:capability','ME virtual-state corrections are not implemented.');
    end
    va=kron(branch.alpha.state,eye(db)); vb=kron(branch.beta.state,eye(db));
    ha=va'*H*va; hb=vb'*H*vb;
    G=-1i*(kron(eye(db),ha)-kron(hb.',eye(db)));
else
    ha=H; hb=H; G=-1i*(kron(eye(d),H)-kron(H.',eye(d)));
end
jump_records=struct('support',{},'rate',{},'left',{},'right',{});

% Embed every wholly included physical jump without changing its rate
for k=1:numel(parameters.jumps)
    jump=parameters.jumps(k);
    if jump.rate==0, continue; end
    included=ismember(jump.support,ids);
    if ~any(included)
        if isfield(context,'exterior_ids')&&any(ismember(jump.support,context.exterior_ids))
            error('cce:exteriorJump','An active exterior jump cannot drive a frozen mean field.');
        end
        continue
    end
    if ~all(included)
        error('cce:jumpSupport','A correlated jump crosses this cluster; group its support.');
    end
    positions=arrayfun(@(id) find(ids==id,1),jump.support);
    rest=setdiff(1:numel(ids),positions,'stable'); perm=[positions rest];
    labels=reshape(1:d,fliplr([1 dims]));
    labels=permute(labels,numel(dims)+2-fliplr([1 perm+1])); labels=labels(:);
    local=kron(jump.operator,eye(prod(dims(rest))));
    op=zeros(d); op(labels,labels)=local;

    % Resolve a closed conditional block or refuse branch leakage
    if conditional
        va=kron(parameters.central.alpha,eye(db));
        vb=kron(parameters.central.beta,eye(db));
        left=va'*op*va; right=vb'*op*vb;
        leakage=norm(op*va-va*left,'fro')+norm(op*vb-vb*right,'fro');

        % Check closure under the anticommutator loss operator as well
        loss=op'*op;
        loss_leakage=norm(loss*va-va*(va'*loss*va),'fro')+...
                     norm(loss*vb-vb*(vb'*loss*vb),'fro');
        if leakage>100*eps*d*norm(op,'fro')||...
           loss_leakage>100*eps*d*norm(loss,'fro')
            error('cce:jumpBranches','A branch-changing jump or loss operator requires megcce.');
        end
    else
        left=op; right=op;
    end

    % Add the native dissipator using tensor-column vectorisation
    dl=size(left,1); a=left'*left; b=right'*right;
    G=G+jump.rate*(kron(conj(right),left)-...
                  (kron(eye(dl),a)+kron(b.',eye(dl)))/2);
    jump_records(end+1)=struct('support',jump.support,'rate',jump.rate,...
                              'left',left,'right',right); %#ok<AGROW>
end

% Return the native-generator convention and exact block identity
meta=struct('physical_ids',ids,'dimension',d,'bath_dimension',db,...
            'central_dimension',d0,'conditional',conditional,...
            'h_alpha',ha,'h_beta',hb,'jumps',jump_records,...
            'vectorization','column','evolution','exp(G*t)',...
            'rate_units','s^-1','hamiltonian',hmeta);

end

% Validate explicit jump supports, ordered matrices and nonnegative rates
function grumble(model,sites,context,parameters)
if ~isstruct(model)||~isfield(model,'physical_ids')||~isstruct(context)
    error('A prepared model and scalar Hamiltonian context are required.');
end
if ~isnumeric(sites)||any(~ismember(sites,model.bath_ids))||...
   numel(unique(sites))~=numel(sites)
    error('sites must contain distinct physical bath IDs.');
end
if ~isstruct(parameters)||~isfield(parameters,'method')||...
   ~ismember(parameters.method,{'mecce','megcce'})||...
   ~isfield(parameters,'jumps')||~isstruct(parameters.jumps)||...
   ~isfield(parameters,'numerics')||~isfield(parameters.numerics,'max_memory_bytes')
    error('Explicit mecce/megcce method, jumps and memory budget are required.');
end
budget=parameters.numerics.max_memory_bytes;
if ~isscalar(budget)||~isreal(budget)||~isfinite(budget)||budget<=0
    error('max_memory_bytes must be positive and finite.');
end
for k=1:numel(parameters.jumps)
    j=parameters.jumps(k);
    if ~all(isfield(j,{'support','operator','rate'}))||...
       ~isnumeric(j.support)||isempty(j.support)||~isrow(j.support)||...
       numel(unique(j.support))~=numel(j.support)||...
       any(~ismember(j.support,[model.central_ids model.bath_ids]))
        error('Every jump requires ordered distinct physical support IDs.');
    end
    n=prod(model.multiplicities(j.support));
    if ~isnumeric(j.operator)||~isequal(size(j.operator),[n n])||...
       any(~isfinite(j.operator),'all')||~isnumeric(j.rate)||...
       ~isscalar(j.rate)||~isreal(j.rate)||~isfinite(j.rate)||j.rate<0
        error('Jump matrices must match support; rates must be finite nonnegative s^-1.');
    end
end
end


