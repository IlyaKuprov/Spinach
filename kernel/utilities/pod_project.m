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
% projected terms and state vectors must remain finite. Projection uses a
% power-of-two scale in column slices bounded by reduced rank before
% matrix products and restores physical units
% afterwards. Unsupported component loss is rejected, not silently removed.
% The Gram tolerance uses a dimension-dependent complex-arithmetic budget
% gamma(8*n)*r, with gamma(k)=k*eps/(1-k*eps).
% Galerkin projection: Benner et al., DOI:10.1137/130932715, section 2.2.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pod_project.m>

function model=pod_project(basis,terms,rho,coil)

% Check Euclidean coordinates and the affine projection contract
grumble(basis,terms,rho,coil);

% Project normalised affine terms and restore their physical units
reduced=cell(size(terms));
for n=1:numel(terms)
    reduced{n}=scaled_projection(basis,terms{n},basis);
end

% Project both state-space vectors with the same range-aware arithmetic
rho=scaled_projection(basis,rho,1); coil=scaled_projection(basis,coil,1);
model=struct('terms',{reduced},'rho',rho,'coil',coil);

end

% Evaluate a Galerkin product in a binary-scaled representation
function projected=scaled_projection(basis,data,right)

% Measure the binary input scale without a second full operator
width=size(basis,2); largest=0;
for n=1:width:size(data,2)
    columns=n:min(n+width-1,size(data,2)); part=data(:,columns);
    largest=max(largest,max(max(abs(real(part)),[],'all'),max(abs(imag(part)),[],'all')));
end
[~,exponent]=log2(largest);

% Allocate only the required projected coordinates while preserving sparse storage
if issparse(basis)&&issparse(data)
    projected=complex(sparse(size(basis,2),size(right,2)));
else
    projected=complex(zeros(size(basis,2),size(right,2)));
end

% Accumulate the Galerkin product from bounded normalised column slices
for n=1:width:size(data,2)
    columns=n:min(n+width-1,size(data,2)); part=data(:,columns);
    scaled=complex(pow2(real(part),-exponent),pow2(imag(part),-exponent));
    if any((real(part)~=0)&(real(scaled)==0),'all')||...
       any((imag(part)~=0)&(imag(scaled)==0),'all')
        error('projection data dynamic range exceeds the normalised double representation.');
    end
    projected=projected+(basis'*scaled)*right(columns,:);
end

% Restore physical units and reject unsupported projected coordinates
projected=complex(pow2(real(projected),exponent),pow2(imag(projected),exponent));
if any(isnan(projected),'all')||any(isinf(projected),'all')
    error('projected terms and states must be finite.');
end

end

% Consistency enforcement
function grumble(basis,terms,rho,coil)
if ~isa(basis,'double')||~ismatrix(basis)||isempty(basis)||...
   any(isnan(basis),'all')||any(isinf(basis),'all')||size(basis,2)>size(basis,1)
    error('basis must be a finite nonempty double n-by-r matrix with r<=n.');
end
roundoff=8*size(basis,1)*eps;
orth_tol=size(basis,2)*roundoff/(1-roundoff);
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


