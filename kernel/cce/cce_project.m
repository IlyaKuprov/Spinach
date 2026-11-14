% Projects a central register onto two explicit states. Syntax:
%
%                   branches=cce_project(model,parameters)
%
% parameters.central.alpha and beta are normalized orthogonal vectors in
% model.central_ids product order. Optional second_order.enabled defaults
% to false; when true, gap_tolerance (rad/s) and overlap_tolerance (0..1)
% are mandatory. Only cce/mecce admit the conditional second-order model.
%
% branches.alpha/beta contain state, projector, energy, expectations and
% susceptibility. For bath operators B_j coupled to central operator S_j,
% the virtual-level correction is sum_jk susceptibility(j,k)*B_j*B_k.
% Operator indices j enumerate x,y,z for successive central spins. No
% eigenvalue-index state selection or bath-denominator approximation is
% hidden: the latter is the stated central-gap perturbation model.
%
% Inputs:
%   model - prepared central-register model.
%   parameters.central.alpha, .beta - normalized orthogonal state vectors.
%   Optional parameters.second_order.enabled; when true,
%   parameters.second_order.gap_tolerance in rad/s and overlap_tolerance are
%   mandatory.
%
% Outputs:
%   branches - alpha/beta state, projector, energy, expectations and
%   susceptibility records, with register operators and projection metadata.
%
% ilya.kuprov@weizmann.ac.il

function branches=cce_project(model,parameters)

% Validate explicit central vectors and perturbation controls
grumble(model,parameters);

% Obtain the isolated central Hamiltonian in the declared physical ordering
[central,meta]=cce_hamiltonian(model,[],struct());
[eigenvectors,eigenvalues]=eig(full(central),'vector');
operators=reshape(meta.operators.',1,[]);

% Resolve the explicitly optional second-order approximation
enabled=isfield(parameters,'second_order')&&parameters.second_order.enabled;
branches.second_order_enabled=enabled; branches.central_hamiltonian=central;
branches.central_ids=model.central_ids; branches.eigenvalues=eigenvalues;
branches.source_hash=model.source_hash;

% Project both user-supplied vectors without selecting arbitrary eigenindices
for name={'alpha','beta'}

    % Record the state and each physical-spin angular momentum expectation
    state=parameters.central.(name{1}); branch.state=state;
    branch.projector=state*state'; branch.energy=real(state'*central*state);
    branch.expectations=zeros(numel(model.central_ids),3);
    for n=1:numel(model.central_ids)
        for axis=1:3
            branch.expectations(n,axis)=real(state'*meta.operators{n,axis}*state);
        end
    end

    % Diagnose overlap and the actual central-state eigen-residual
    overlaps=abs(eigenvectors'*state).^2;
    [branch.state_overlap,selected]=max(overlaps);
    branch.eigen_residual=norm(central*state-branch.energy*state);
    branch.minimum_gap=min(abs(eigenvalues((1:numel(eigenvalues))~=selected)-...
                               eigenvalues(selected)));
    susceptibility=zeros(numel(operators),numel(operators));

    % Sum virtual central transitions only when their denominators are resolved
    if enabled
        tolerance=parameters.second_order;
        if branch.state_overlap<1-tolerance.overlap_tolerance
            error('central.%s is not an isolated central eigenstate within overlap_tolerance.',name{1});
        end
        if branch.minimum_gap<=tolerance.gap_tolerance
            error('central.%s has a near-degenerate second-order denominator.',name{1});
        end
        for n=1:numel(eigenvalues)
            if n==selected, continue; end
            matrix_elements=zeros(numel(operators),1);
            for k=1:numel(operators)
                matrix_elements(k)=state'*operators{k}*eigenvectors(:,n);
            end
            susceptibility=susceptibility+...
                (matrix_elements*matrix_elements')/(eigenvalues(selected)-eigenvalues(n));
        end
    end

    % Store branch-specific first- and second-order coefficients
    branch.susceptibility=susceptibility; branches.(name{1})=branch;
end

end

% Reject invalid state vectors and unresolved perturbation settings
function grumble(model,parameters)
if ~isstruct(model)||~all(isfield(model,{'central_ids','multiplicities','source_hash'}))
    error('model must be returned by cce_prepare.');
end
if ~isstruct(parameters)||~isfield(parameters,'central')||...
   ~isstruct(parameters.central)||~all(isfield(parameters.central,{'alpha','beta'}))
    error('central.alpha and central.beta must be explicitly supplied.');
end
dimension=prod(model.multiplicities(model.central_ids));
for name={'alpha','beta'}
    state=parameters.central.(name{1});
    if ~isa(state,'double')||~isequal(size(state),[dimension 1])||any(~isfinite(state))
        error('central.%s must be a finite double column of central dimension.',name{1});
    end
    if abs(norm(state)-1)>dimension*eps
        error('central.%s must have unit norm.',name{1});
    end
end
if abs(parameters.central.alpha'*parameters.central.beta)>dimension*eps
    error('central.alpha and central.beta must be orthogonal.');
end
if isfield(parameters,'second_order')
    control=parameters.second_order;
    if ~isstruct(control)||~isscalar(control)||~isfield(control,'enabled')||...
       ~islogical(control.enabled)||~isscalar(control.enabled)||...
       ~all(ismember(fieldnames(control),{'enabled','gap_tolerance','overlap_tolerance'}))
        error('second_order requires a scalar logical enabled field.');
    end
    if control.enabled
        if ~isfield(parameters,'method')||~ismember(parameters.method,{'cce','mecce'})
            error('second-order projection is defined only for cce and mecce.');
        end
        if ~all(isfield(control,{'gap_tolerance','overlap_tolerance'}))
            error('second-order gap and overlap tolerances must be explicit.');
        end
        for name={'gap_tolerance','overlap_tolerance'}
            value=control.(name{1});
            if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<0
                error('invalid second_order.%s.',name{1});
            end
        end
        if control.gap_tolerance==0||control.overlap_tolerance>=1
            error('positive gap tolerance and overlap tolerance below one required.');
        end
    elseif isfield(control,'gap_tolerance')||isfield(control,'overlap_tolerance')
        error('second-order tolerances must be absent when correction is disabled.');
    end
end
end

