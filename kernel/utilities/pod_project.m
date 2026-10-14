% Fixed-basis Galerkin projection of affine Spinach Liouvillian terms,
% an initial state, and a detection vector. Syntax:
%
%                     model=pod_project(basis,terms,rho,coil)
%
% Parameters:
%
%     basis - finite double n-by-r orthonormal basis in Euclidean
%             Liouville state coordinates, with r<=n
%
%     terms - nonempty cell column of finite double n-by-n Liouvillian
%             terms, sparse or dense, all in the same state coordinates
%
%     rho   - finite double n-by-1 initial state
%
%     coil  - finite double n-by-1 detection vector; signal is coil'*rho
%
% Outputs:
%
%     model.terms - projected affine terms, basis'*terms{k}*basis
%     model.rho   - projected initial state, basis'*rho
%     model.coil  - projected detection vector, basis'*coil
%
% Use L=H+1i*R+1i*K in rad/s with time in seconds and exp(-1i*L*t).
% Projection commutes algebraically with affine assembly, but does not
% guarantee accurate reduced trajectories. Terms, states, and basis must
% share ordering and projectors; no coordinate conversion is inferred.
% Non-Euclidean metrics require explicit caller-supplied whitening first.
% Orthogonality alone does not preserve positivity or trace. Finiteness
% checks preserve sparse zeros without extracting dense entries. Sparse
% bases retain a sparse identity in the orthogonality check. Computed
% projected terms and state vectors must remain finite.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pod_project.m>

function model=pod_project(basis,terms,rho,coil)

% Check the projection contract
grumble(basis,terms,rho,coil);

% Project the fixed affine terms and both state-space vectors
reduced=cell(size(terms));
for n=1:numel(terms)
    reduced{n}=basis'*terms{n}*basis;
    if any(isnan(reduced{n}),'all')||any(isinf(reduced{n}),'all')
        error('projected Liouvillian terms must be finite.');
    end
end
rho=basis'*rho; coil=basis'*coil;
if any(isnan(rho),'all')||any(isinf(rho),'all')||...
   any(isnan(coil),'all')||any(isinf(coil),'all')
    error('projected states must be finite.');
end
model=struct('terms',{reduced},'rho',rho,'coil',coil);

end

% Consistency enforcement
function grumble(basis,terms,rho,coil)
if ~isa(basis,'double')||~ismatrix(basis)||isempty(basis)||...
   any(isnan(basis),'all')||any(isinf(basis),'all')||size(basis,2)>size(basis,1)
    error('basis must be a finite nonempty double n-by-r matrix with r<=n.');
end
orth_tol=size(basis,1)*size(basis,2)*eps;
if norm(basis'*basis-speye(size(basis,2)),'fro')>orth_tol
    error('basis must be orthonormal in Euclidean state coordinates.');
end
if ~iscell(terms)||~iscolumn(terms)||isempty(terms)
    error('terms must be a nonempty cell column.');
end
for n=1:numel(terms)
    if ~isa(terms{n},'double')||~isequal(size(terms{n}),[size(basis,1) size(basis,1)])||...
       any(isnan(terms{n}),'all')||any(isinf(terms{n}),'all')
        error('each term must be a finite double n-by-n matrix.');
    end
end
if ~isa(rho,'double')||~isequal(size(rho),[size(basis,1) 1])||any(isnan(rho),'all')||any(isinf(rho),'all')
    error('rho must be a finite double n-by-1 column.');
end
if ~isa(coil,'double')||~isequal(size(coil),[size(basis,1) 1])||any(isnan(coil),'all')||any(isinf(coil),'all')
    error('coil must be a finite double n-by-1 column.');
end
end


