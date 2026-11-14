% Extract central density elements without discarding complex phase. Syntax:
%
%                  observed=cce_observe(rho,parameters)
%
% rho is d-by-d-by-ntime. quantity is density, coherence, or observable.
% Coherence uses central.alpha/beta; observable uses operators (a matrix
% or cell of matrices). normalized is an explicit logical flag for scalar
% observables; its reference is central.rho, never the first output time.
% Density elements are flattened in MATLAB column-major order.
%
% Inputs:
%   rho - central density array, d-by-d-by-ntime.
%   parameters - quantity=density|coherence|observable; central alpha/beta
%   vectors or observable operators (matrix/cell); explicit normalized for
%   scalar quantities.
%
% Outputs:
%   observed - numeric time-by-observable complex array; unresolved scalar
%   normalization raises cce_observe:initial_zero.
%
% See also cce_cluster, cce_master.
%
% ilya.kuprov@weizmann.ac.il

function observed=cce_observe(rho,parameters)

% Validate the density stack and requested observable
grumble(rho,parameters);

% Extract each time point in an explicitly documented element order
d=size(rho,1); nt=size(rho,3);
switch parameters.quantity
    case 'density'
        observed=reshape(rho,d*d,nt).';
    case 'coherence'
        observed=complex(zeros(nt,1));
        for k=1:nt
            observed(k)=parameters.central.alpha'*rho(:,:,k)*parameters.central.beta;
        end
    case 'observable'
        operators=parameters.operators;
        if ~iscell(operators), operators={operators}; end
        observed=complex(zeros(nt,numel(operators)));
        for k=1:nt
            for j=1:numel(operators)
                observed(k,j)=trace(operators{j}*rho(:,:,k));
            end
        end
end

% Normalize only against the explicitly supplied initial central state
if isfield(parameters,'normalized')&&parameters.normalized
    if strcmp(parameters.quantity,'density')
        error('cce_observe:normalization','A density matrix is not normalized elementwise.');
    end
    reference_parameters=parameters; reference_parameters.normalized=false;
    reference=cce_observe(parameters.central.rho,reference_parameters);
    bound=size(rho,1)^2*eps*max(1,norm(parameters.central.rho,'fro'));
    % Propagate density uncertainty in the units of each observable
    if strcmp(parameters.quantity,'observable')
        bound=bound*cellfun(@(op) norm(op,'fro'),operators);
    end
    if any(abs(reference)<=bound)
        error('cce_observe:initial_zero','Initial observable is unresolved or zero.');
    end
    observed=observed./reference;
end

end

% Validate dimensions, finite values and observation definitions
function grumble(rho,parameters)
if ~isnumeric(rho)||ndims(rho)>3||size(rho,1)~=size(rho,2)||...
   isempty(rho)||any(~isfinite(rho),'all')
    error('rho must be a finite square matrix or stack.');
end
if ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'quantity')||...
   ~ismember(parameters.quantity,{'density','coherence','observable'})
    error('quantity must be density, coherence, or observable.');
end
if strcmp(parameters.quantity,'coherence')
    if ~isfield(parameters,'central')||...
       ~all(isfield(parameters.central,{'alpha','beta'}))
        error('Coherence requires central.alpha and central.beta.');
    end
    for name={'alpha','beta'}
        state=parameters.central.(name{1});
        if ~isnumeric(state)||~isequal(size(state),[size(rho,1) 1])||...
           any(~isfinite(state))||abs(norm(state)-1)>size(rho,1)*eps
            error('Central observation states must be normalized columns.');
        end
    end
    if abs(parameters.central.alpha'*parameters.central.beta)>size(rho,1)*eps
        error('Central observation states must be orthogonal.');
    end
end
if strcmp(parameters.quantity,'observable')
    if ~isfield(parameters,'operators'), error('operators must be supplied.'); end
    operators=parameters.operators;
    if ~iscell(operators), operators={operators}; end
    if isempty(operators), error('operators must be nonempty.'); end
    for k=1:numel(operators)
        if ~isnumeric(operators{k})||~isequal(size(operators{k}),size(rho,[1 2]))||...
           any(~isfinite(operators{k}),'all')
            error('Each operator must match the central dimension.');
        end
    end
end
if isfield(parameters,'normalized')
    if ~islogical(parameters.normalized)||~isscalar(parameters.normalized)
        error('normalized must be a logical scalar.');
    end
    if parameters.normalized&&(~isfield(parameters,'central')||...
       ~isfield(parameters.central,'rho'))
        error('Normalization requires central.rho.');
    end
end
end


