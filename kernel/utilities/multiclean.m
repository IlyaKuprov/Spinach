% Memory-prioritised matrix multiplication with immediate clean_up of
% each output panel. No inner-product partial sum is rounded. Syntax:
%
%                         C=multiclean(A,B,tol)
%
% Parameters:
%
%     A,B  - two-dimensional single or double matrices, full or sparse,
%            on CPU or GPU; dimensions must conform for A*B, including
%            MATLAB scalar multiplication
%
%     tol  - real, finite, non-negative single or double CPU scalar;
%            zero returns native A*B without panel multiplication
%
% Outputs:
%
%     C    - cleaned product, with the precision of clean_up(A*B) and
%            the device of A*B (single tol can lower double precision)
%
% The three-argument interface preserves the native A*B sparse/full
% representation; the caller retains its spin-system storage policy.
% Rounding is the existing tol*round((1/tol)*panel) operation, not
% magnitude thresholding. With tol=0, A*B is returned directly.
% Dense and sparse inputs use one multiply path; partial inner products
% are never rounded.
%
% Panels have at most one input-storage footprint, counting three product
% buffers and the sliced right operand conservatively as dense. On GPUs,
% this budget is also limited by current AvailableMemory. This bounds the
% explicit uncleaned intermediates, not undocumented BLAS/CUDA workspace
% or the cleaned output. Sparse assembly temporarily holds triplets as
% well as the result. A dense result must itself fit in available memory.
% GPU sparse inputs require MATLAB R2026b contiguous-submatrix indexing.
% Blocked BLAS reductions can differ in their last bits from A*B, hence
% values on a rounding boundary need not give bitwise-identical results.
%
% ilya.kuprov@weizmann.ac.il

function C=multiclean(A,B,tol)

% Check consistency
grumble(A,B,tol);

% Disabled clean-up must also preserve native multiplication arithmetic
if tol==0, C=A*B; return; end

% Obtain product metadata without computing an uncleaned product
if isscalar(A)
    proto=A*zeros(0,0,'like',B); nrows=size(B,1); ncols=size(B,2);
elseif isscalar(B)
    proto=zeros(0,0,'like',A)*B; nrows=size(A,1); ncols=size(A,2);
else
    proto=zeros(0,0,'like',A)*zeros(0,0,'like',B);
    nrows=size(A,1); ncols=size(B,2);
end

% Specify only the local policy needed by clean_up
spin_system.sys.disable={};
if issparse(proto)
    spin_system.tols.dense_matrix=Inf;
else
    spin_system.tols.dense_matrix=0;
end
spin_system.tols.small_matrix=0;

% Include clean_up scalar arithmetic in the output precision metadata
clean_proto=clean_up(spin_system,proto,tol);

% Empty products require no arithmetic or storage conversion
if (nrows==0)||(ncols==0)
    C=zeros(nrows,ncols,'like',clean_proto); return
end

% Preserve matrix arithmetic by placing a scalar on the left
if isscalar(B)&&~isscalar(A)
    C=multiclean(B,A,tol); return
end

% Measure input storage, not the size of the hypothetical dense product
if isa(proto,'gpuArray')
    input_bytes=0;
    for input={A,B}
        data=input{1};
        if isa(data,'gpuArray'), precision=classUnderlying(data);
        else, precision=class(data); end
        data_bytes=8-4*strcmp(precision,'single');
        data_bytes=data_bytes*(1+~isreal(data));
        if issparse(data)
            input_bytes=input_bytes+(data_bytes+8)*double(gather(nnz(data)))+...
                        8*(size(data,2)+1);
        else
            input_bytes=input_bytes+data_bytes*numel(data);
        end
    end
    device=gpuDevice; input_bytes=min(input_bytes,device.AvailableMemory);
    precision=classUnderlying(proto);
else
    info=whos('A','B'); input_bytes=sum([info.bytes]);
    precision=class(proto);
end

% Bound live panel storage independently of output sparsity
value_bytes=(8-4*strcmp(precision,'single'))*(1+(~isreal(A)||~isreal(B)));
entry_bytes=value_bytes+8;
slice_bytes=8*size(B,1)*(1+~isreal(B)+issparse(B))+8*issparse(B);
col_bytes=3*(nrows*entry_bytes+8)+slice_bytes;
panel_width=max(1,min(ncols,floor(input_bytes/max(1,col_bytes))));
if ncols>1, panel_width=min(panel_width,ceil(ncols/2)); end

% Retain only cleaned entries, without any caller-specific density policy
panel_count=ceil(ncols/panel_width);
dense_output=~issparse(proto);
if dense_output
    C=zeros(nrows,ncols,'like',clean_proto);
else
    row_parts=cell(1,panel_count); col_parts=cell(1,panel_count);
    val_parts=cell(1,panel_count);
end

% Multiply complete inner products and immediately clean each panel
for panel_index=1:panel_count
    first_col=1+(panel_index-1)*panel_width;
    last_col=min(ncols,first_col+panel_width-1);
    panel=A*B(:,first_col:last_col);
    panel=clean_up(spin_system,panel,tol);

    % Sparse GPU arrays cannot be assigned by subscript
    if dense_output
        C(:,first_col:last_col)=full(panel);
    else
        [row_index,col_index,values]=find(panel);
        row_parts{panel_index}=row_index(:);
        col_parts{panel_index}=col_index(:);
        val_parts{panel_index}=values(:);
        col_parts{panel_index}=col_parts{panel_index}+first_col-1;
    end
    clear panel
end

% Sparse GPU arrays cannot be concatenated; assemble retained triplets
if ~dense_output
    row_index=vertcat(row_parts{:}); clear row_parts
    col_index=vertcat(col_parts{:}); clear col_parts
    values=vertcat(val_parts{:}); clear val_parts
    C=sparse(row_index,col_index,values,nrows,ncols);
end

end

% Consistency enforcement
function grumble(A,B,tol)
if (~isnumeric(A))||(~isfloat(A))||(~ismatrix(A))
    error('A must be a two-dimensional single or double array.');
end
if (~isnumeric(B))||(~isfloat(B))||(~ismatrix(B))
    error('B must be a two-dimensional single or double array.');
end
if (~isscalar(A))&&(~isscalar(B))&&(size(A,2)~=size(B,1))
    error('matrix dimensions must conform for A*B.');
end
if isa(tol,'gpuArray')||(~isnumeric(tol))||(~isfloat(tol))||(~isreal(tol))||...
   (~isscalar(tol))||(~isfinite(tol))||(tol<0)
    error('tol must be a finite non-negative real single or double CPU scalar.');
end
end


