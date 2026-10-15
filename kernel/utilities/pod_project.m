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
% afterwards. The fixed Hermitian basis and its component range are prepared
% once, outside the projection slices. Dense range scans use bounded columns.
% Real operands retain real slice and output storage; complex products promote
% storage only when required. Unsupported component loss is rejected, not removed.
% Nonzero operand-component pairs must remain nonzero at the working scale;
% this conservative product-range contract can reject a finite physical result.
% The Gram tolerance uses a dimension-dependent complex-arithmetic budget
% gamma(8*n)*r, with gamma(k)=k*eps/(1-k*eps).
% Galerkin projection: Benner et al., DOI:10.1137/130932715, section 2.2.
%
% Physical restoration rejects a nonzero working component rounded to zero,
% including loss of just one component of a complex projected coordinate.
% Representable subnormal coordinates remain supported.
%
% ilya.kuprov@weizmann.ac.il
% arnab.acharya@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pod_project.m>

function model=pod_project(basis,terms,rho,coil)

% Check Euclidean coordinates and the affine projection contract
grumble(basis,terms,rho,coil);

% Prepare the fixed Hermitian basis and its component range once
basis_min=component_min(basis); left_basis=basis';

% Project normalised affine terms and restore their physical units
reduced=cell(size(terms));
for n=1:numel(terms)
    reduced{n}=scaled_projection(left_basis,basis_min,terms{n},basis);
end

% Project both state-space vectors with the same range-aware arithmetic
rho=scaled_projection(left_basis,basis_min,rho,1);
coil=scaled_projection(left_basis,basis_min,coil,1);
model=struct('terms',{reduced},'rho',rho,'coil',coil);

end

% Evaluate a Galerkin product in a binary-scaled representation
function projected=scaled_projection(left_basis,basis_min,data,right)

% Measure the binary input scale without a second full operator
width=size(left_basis,1); largest=0;
for n=1:width:size(data,2)
    columns=n:min(n+width-1,size(data,2)); part=data(:,columns);
    largest=max(largest,max(max(abs(real(part)),[],'all'),max(abs(imag(part)),[],'all')));
end
[~,exponent]=log2(largest);

% Allocate only the required projected coordinates while preserving sparse storage
if issparse(left_basis)&&issparse(data)
    projected=sparse(size(left_basis,1),size(right,2));
else
    projected=zeros(size(left_basis,1),size(right,2));
end

% Accumulate the Galerkin product from bounded normalised column slices
for n=1:width:size(data,2)
    columns=n:min(n+width-1,size(data,2)); part=data(:,columns);
    scaled=pow2(real(part),-exponent);
    if ~isreal(part), scaled=complex(scaled,pow2(imag(part),-exponent)); end
    if any((real(part)~=0)&(real(scaled)==0),'all')||...
       (~isreal(part)&&any((imag(part)~=0)&(imag(scaled)==0),'all'))
        error('projection data dynamic range exceeds the normalised double representation.');
    end
    if basis_min*component_min(scaled)==0
        error('projection product dynamic range exceeds the normalised double representation.');
    end
    projected=projected+checked_product(left_basis*scaled,right(columns,:));
end

% Restore physical units while tracking every nonzero real and imaginary component
real_count=nnz(real(projected));
if isreal(projected)
    imag_count=0; projected=pow2(projected,exponent);
else
    imag_count=nnz(imag(projected));
    projected=complex(pow2(real(projected),exponent),pow2(imag(projected),exponent));
end
if any(isnan(projected),'all')||any(isinf(projected),'all')
    error('projected terms and states must be finite.');
end
if nnz(real(projected))~=real_count||...
   (imag_count~=0&&nnz(imag(projected))~=imag_count)
    error('projected coordinate dynamic range exceeds the physical double representation.');
end

end

% Reject unsupported operand ranges before a bounded matrix product
function product=checked_product(left,right)

% Bound every nonzero real component product without forming a product array
if component_min(left)*component_min(right)==0
    error('projection product dynamic range exceeds the normalised double representation.');
end
product=left*right;

end

% Find the smallest nonzero component with bounded dense column storage
function minimum=component_min(data)

% Preserve stored sparse traversal without copying dense matrix components
minimum=Inf;
if issparse(data)
    parts=nonzeros(real(data));
    if ~isempty(parts), minimum=min(abs(parts)); end
    parts=nonzeros(imag(data));
    if ~isempty(parts), minimum=min(minimum,min(abs(parts))); end
else
    for k=1:size(data,2)
        parts=real(data(:,k)); parts=parts(parts~=0);
        if ~isempty(parts), minimum=min(minimum,min(abs(parts))); end
        parts=imag(data(:,k)); parts=parts(parts~=0);
        if ~isempty(parts), minimum=min(minimum,min(abs(parts))); end
    end
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


