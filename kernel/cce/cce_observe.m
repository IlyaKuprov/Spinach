% Extract observed density elements without discarding complex phase. Syntax:
%
%                  observed=cce_observe(rho,parameters)
%
% rho is d-by-d-by-ntime. quantity is density, signal, or observable.
% Signal uses transition.alpha/beta; observable uses coil (a matrix
% or cell of matrices), with signal trace(coil'*rho). normalized is an
% explicit logical flag for scalar
% observables; its reference is rho0, never the first output time.
% Density elements are flattened in MATLAB column-major order.
%
% Inputs:
%   rho - observed density array, d-by-d-by-ntime.
%   parameters - quantity=density|signal|observable; transition.alpha/beta
%   vectors or receiver coil (matrix/cell); explicit normalized for
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

% Reject unsupported optional inputs before selecting an execution path
cce_validate(parameters);

% Validate the density stack and requested observable
grumble(rho,parameters);

% Extract each time point in an explicitly documented element order
d=size(rho,1); nt=size(rho,3);
switch parameters.quantity
    case 'density'
        observed=reshape(rho,d*d,nt).';
    case 'signal'
        observed=complex(zeros(nt,1));
        for k=1:nt
            observed(k)=parameters.transition.alpha'*rho(:,:,k)*parameters.transition.beta;
        end
    case 'observable'
        coils=parameters.coil;
        if ~iscell(coils), coils={coils}; end
        observed=complex(zeros(nt,numel(coils)));
        for k=1:nt
            for j=1:numel(coils)
                observed(k,j)=trace(coils{j}'*rho(:,:,k));
            end
        end
end

% Normalize only against the explicitly supplied initial observed state
if isfield(parameters,'normalized')&&parameters.normalized
    if strcmp(parameters.quantity,'density')
        error('cce_observe:normalization','A density matrix is not normalized elementwise.');
    end
    reference_parameters=parameters; reference_parameters.normalized=false;
    reference=cce_observe(parameters.rho0,reference_parameters);
    bound=size(rho,1)^2*eps*max(1,norm(parameters.rho0,'fro'));
    % Propagate density uncertainty in the units of each observable
    if strcmp(parameters.quantity,'observable')
        bound=bound*cellfun(@(op) norm(op,'fro'),coils(:).');
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
   ~ismember(parameters.quantity,{'density','signal','observable'})
    error('quantity must be density, signal, or observable.');
end
if strcmp(parameters.quantity,'signal')
    if ~isfield(parameters,'transition')||...
       ~all(isfield(parameters.transition,{'alpha','beta'}))
        error('Signal detection requires parameters.transition.alpha and parameters.transition.beta.');
    end
    for name={'alpha','beta'}
        state=parameters.transition.(name{1});
        if ~isnumeric(state)||~isequal(size(state),[size(rho,1) 1])||...
           any(~isfinite(state))||abs(norm(state)-1)>size(rho,1)*eps
            error('Observed-spin observation states must be normalized columns.');
        end
    end
    if abs(parameters.transition.alpha'*parameters.transition.beta)>size(rho,1)*eps
        error('Observed-spin observation states must be orthogonal.');
    end
end
if strcmp(parameters.quantity,'observable')
    if ~isfield(parameters,'coil'), error('coil must be supplied.'); end
    coils=parameters.coil;
    if ~iscell(coils), coils={coils}; end
    if isempty(coils), error('coils must be nonempty.'); end
    for k=1:numel(coils)
        if ~isnumeric(coils{k})||~isequal(size(coils{k}),size(rho,[1 2]))||...
           any(~isfinite(coils{k}),'all')
            error('Each receiver coil must match the observed-spin dimension.');
        end
    end
end
if isfield(parameters,'normalized')
    if ~islogical(parameters.normalized)||~isscalar(parameters.normalized)
        error('normalized must be a logical scalar.');
    end
    if parameters.normalized&&~isfield(parameters,'rho0')
        error('Normalization requires rho0.');
    end
end
end


