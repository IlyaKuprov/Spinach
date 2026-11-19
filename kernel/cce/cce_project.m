% Projects an observed-spin system onto two explicit states. Syntax:
%
%                   manifolds=cce_project(model,parameters)
%
% parameters.transition.alpha and beta are normalised orthogonal vectors in
% model.observed_ids product order. Optional second_order.enabled defaults
% to false; when true, gap_tolerance (rad/s) and overlap_tolerance (0..1)
% are mandatory. Only cce/mecce admit the projected second-order model.
%
% manifolds.alpha/beta contain vector, projector, energy, expectations, and
% susceptibility. For bath operators B_j coupled to observed operator S_j,
% the virtual-level correction is sum_jk susceptibility(j,k)*B_j*B_k.
% Operator indices j enumerate x,y,z for successive observed spins. No
% eigenvalue-index state selection or bath-denominator approximation is
% hidden: the latter is the stated observed-gap perturbation model.
%
% Parameters:
%
%   model - prepared observed-spin system model.
%   parameters.transition.alpha, .beta - normalised orthogonal state vectors.
%   Optional parameters.second_order.enabled; when true,
%   parameters.second_order.gap_tolerance in rad/s and overlap_tolerance are
%   mandatory.
%
% Outputs:
%
%   manifolds - alpha/beta state, projector, energy, expectations, and
%   susceptibility records, with observed-spin system operators and projection metadata.
%
% See also cce_compile, cce_generators.
%
% ilya.kuprov@weizmann.ac.il

function manifolds=cce_project(model,parameters)

% Validate explicit observed vectors and perturbation controls
grumble(model,parameters);

% Obtain the isolated observed Hamiltonian in the declared physical ordering
[H_obs,meta]=cce_hamiltonian(model,[],struct());
[eigenvectors,eigenvalues]=eig(full(H_obs),'vector');
operators=reshape(meta.operators.',1,[]);

% Resolve the explicitly optional second-order approximation
enabled=isfield(parameters,'second_order')&&parameters.second_order.enabled;
manifolds.second_order_enabled=enabled;
manifolds.observed_hamiltonian=H_obs;
manifolds.observed_ids=model.observed_ids;
manifolds.eigenvalues=eigenvalues;
manifolds.source_hash=model.source_hash;

% Project both user-supplied vectors without selecting arbitrary eigenindices
for name={'alpha','beta'}

    % Record the state and each physical-spin angular momentum expectation
    state=parameters.transition.(name{1});
    manifold.vector=state;
    manifold.projector=state*state';
    manifold.energy=real(state'*H_obs*state);
    manifold.expectations=zeros(numel(model.observed_ids),3);
    for n=1:numel(model.observed_ids)
        for a=1:3
            manifold.expectations(n,a)=real(state'*meta.operators{n,a}*state);
        end
    end

    % Diagnose overlap and the actual observed-state eigen-residual
    overlaps=abs(eigenvectors'*state).^2;
    [manifold.state_overlap,selected]=max(overlaps);
    manifold.eigen_residual=norm(H_obs*state-manifold.energy*state);
    manifold.minimum_gap=min(abs(eigenvalues((1:numel(eigenvalues))~=selected)-...
                               eigenvalues(selected)));
    susceptibility=zeros(numel(operators),numel(operators));

    % Sum virtual observed transitions only when their denominators are resolved
    if enabled
        tolerance=parameters.second_order;
        if manifold.state_overlap<1-tolerance.overlap_tolerance
            error('transition.%s is not an isolated observed eigenstate within overlap_tolerance.',name{1});
        end
        if manifold.minimum_gap<=tolerance.gap_tolerance
            error('transition.%s has a near-degenerate second-order denominator.',name{1});
        end
        for n=1:numel(eigenvalues)
            if n==selected
                continue;
            end
            matrix_elements=zeros(numel(operators),1);
            for k=1:numel(operators)
                matrix_elements(k)=state'*operators{k}*eigenvectors(:,n);
            end
            susceptibility=susceptibility+...
                (matrix_elements*matrix_elements')/(eigenvalues(selected)-eigenvalues(n));
        end
    end

    % Store manifold-specific first- and second-order coefficients
    manifold.susceptibility=susceptibility;
    manifolds.(name{1})=manifold;

end

end

% Reject invalid state vectors and unresolved perturbation settings
function grumble(model,parameters)
if ~isstruct(model)||~all(isfield(model,{'observed_ids','multiplicities','source_hash'}))
    error('model must be returned by cce_compile.');
end
if ~isstruct(parameters)||~isfield(parameters,'transition')||...
   ~isstruct(parameters.transition)||~all(isfield(parameters.transition,{'alpha','beta'}))
    error('transition.alpha and transition.beta must be explicitly supplied.');
end
dimension=prod(model.multiplicities(model.observed_ids));
for name={'alpha','beta'}
    state=parameters.transition.(name{1});
    if ~isa(state,'double')||~isequal(size(state),[dimension 1])||any(~isfinite(state))
        error('transition.%s must be a finite double column of observed dimension.',name{1});
    end
    if abs(norm(state)-1)>dimension*eps
        error('transition.%s must have unit norm.',name{1});
    end
end
if abs(parameters.transition.alpha'*parameters.transition.beta)>dimension*eps
    error('transition.alpha and transition.beta must be orthogonal.');
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


