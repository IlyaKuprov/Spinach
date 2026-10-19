% Memory-prioritised matrix multiplication with panel-wise cleanup. Syntax:
%
%                 C=multiclean(spin_system,A,B,tol)
%
% A,B are conforming matrices on CPU or GPU; tol is a finite, real,
% non-negative floating-point scalar. Zero leaves entries unrounded.
% Storage and rounding decisions belong to clean_up. Sparse GPU products
% use cuda_sparse_by_sparse internally. Panel assembly may change floating-
% point roundoff relative to a single native multiplication.
%
% ilya.kuprov@weizmann.ac.il

function C=multiclean(spin_system,A,B,tol)

% Check consistency
grumble(A,B,tol);

% Scalars, empty products and symbolic objects do not need panel storage
if isscalar(A)||isscalar(B)||isempty(A)||isempty(B)||...
   (isobject(A)&&~isa(A,'gpuArray'))||(isobject(B)&&~isa(B,'gpuArray'))
    C=clean_up(spin_system,A*B,tol); return
end

% Put both operands on the execution device
if isa(A,'gpuArray')||isa(B,'gpuArray')
    A=gpuArray(A); B=gpuArray(B);
end

% Budget three dense panel arrays against the input nonzero count
width=max(1,floor((nnz(A)+nnz(B))/(3*size(A,1))));
width=min(width,size(B,2));
nblocks=ceil(size(B,2)/width); sparse_product=issparse(A)&&issparse(B);
if sparse_product, rows=cell(1,nblocks); cols=rows; vals=rows; end

% Multiply complete columns, clean them and retain their sparse storage
for n=1:nblocks
    columns=(1+(n-1)*width):min(n*width,size(B,2));
    if isa(A,'gpuArray')&&issparse(A)&&issparse(B)&&...
       strcmp(classUnderlying(A),'double')&&strcmp(classUnderlying(B),'double')
        panel=cuda_sparse_by_sparse(A,B(:,columns));
    else
        panel=A*B(:,columns);
    end
    panel=clean_up(spin_system,panel,tol);
    if sparse_product
        [rows{n},cols{n},vals{n}]=find(panel);
        rows{n}=rows{n}(:); cols{n}=cols{n}(:)+columns(1)-1; vals{n}=vals{n}(:);
    else
        if n==1, C=zeros(size(A,1),size(B,2),'like',full(panel)); end
        C(:,columns)=full(panel);
    end
end

% Let clean_up choose storage for the assembled matrix
if sparse_product
    C=sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),size(A,1),size(B,2));
    clear rows cols vals
end
C=clean_up(spin_system,C,tol);

end

% Consistency enforcement
function grumble(A,B,tol)
if (~ismatrix(A))||(~ismatrix(B))
    error('A and B must be matrices.');
end
if (~isscalar(A))&&(~isscalar(B))&&(size(A,2)~=size(B,1))
    error('matrix dimensions must conform for A*B.');
end
if (~isfloat(tol))||(~isreal(tol))||(~isscalar(tol))||...
   (~isfinite(tol))||(tol<0)
    error('tol must be a finite non-negative real floating-point scalar.');
end
end


