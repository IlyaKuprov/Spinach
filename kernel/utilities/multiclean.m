% Matrix multiplication with clean-up. Syntax:
%
%                 C=multiclean(spin_system,A,B,tol)
%
% Parameters:
%
%     spin_system - Spinach system structure supplying cleanup
%                   policy and the parallel worker count
%
%     A,B         - conforming matrices on CPU or GPU; scalar and
%                   symbolic-operator products are also accepted
%
%     tol         - finite, real, non-negative floating-point scalar;
%                   zero leaves entries unrounded
%
% Outputs:
%
%     C           - cleaned product with the dimensions of A*B
%
% When the full product does not fit into the memory share of the
% current worker (host free memory divided by the worker count), it
% is computed in column panels that are cleaned up as they are made,
% which limits the fill-in of sparse products. GPU products always
% use ordinary multiplication; sparse double GPU products use the
% cuda_sparse_by_sparse backend. Panel multiplication may differ
% from one native product by numerical roundoff.
%
% ilya.kuprov@weizmann.ac.il

function C=multiclean(spin_system,A,B,tol)

% Check consistency
grumble(A,B,tol);

% GPU arrays and symbolic objects: ordinary cleaned-up multiplication
if isobject(A)||isobject(B)
    if isa(A,'gpuArray')&&isa(B,'gpuArray')&&issparse(A)&&issparse(B)&&...
       strcmp(classUnderlying(A),'double')&&strcmp(classUnderlying(B),'double')
        C=cuda_sparse_by_sparse(A,B);
    else
        C=A*B;
    end
    C=clean_up(spin_system,C,tol); return
end

% Memory share of this worker, assuming all workers do the same job
budget=free_memory()/max(1,spin_system.sys.parallel{2});

% Bytes per full product element, with two clean_up temporaries
elem=full(ones(1,'like',A)*ones(1,'like',B)); elem=whos('elem');
elem_bytes=3*elem.bytes*(1+(~isreal(A)||~isreal(B)));

% Full product fits into the memory share: ordinary multiplication
if size(A,1)*size(B,2)*elem_bytes<=budget
    C=clean_up(spin_system,A*B,tol); return
end

% Column panels of the product that fit into the memory share
width=max(1,floor(budget/(size(A,1)*elem_bytes)));
C=cell(1,ceil(size(B,2)/width));
for n=1:numel(C)
    columns=(1+(n-1)*width):min(n*width,size(B,2));
    C{n}=sparse(clean_up(spin_system,A*B(:,columns),tol));
end
C=clean_up(spin_system,[C{:}],tol);

end

% Available host memory in bytes
function bytes=free_memory()
if ispc
    stats=memory(); bytes=stats.MemAvailableAllArrays;
elseif ismac
    [~,txt]=system('vm_stat'); pages=regexp(txt,'Pages (free|inactive):\s+(\d+)','tokens');
    bytes=4096*sum(cellfun(@(c)str2double(c{2}),pages));
else
    txt=fileread('/proc/meminfo');
    bytes=1024*str2double(regexp(txt,'MemAvailable:\s+(\d+)','tokens','once'));
end
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


