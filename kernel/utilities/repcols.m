% Replicates selected columns of a sparse matrix. Syntax:
%
%                     B=repcols(A,column_numbers,replication_counts)
%
% Parameters:
%
%   A                  - m-by-n sparse matrix
%   column_numbers     - vector of column indices to replicate
%   replication_counts - vector of positive integers specifying how many
%                        copies of each column to make
%
% Output:
%
%   B                  - m-by-(n+sum(replication_counts-1)) sparse matrix
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=repcols.m>

function B=repcols(A,column_numbers,replication_counts)

% Check consistency
grumble(A,column_numbers,replication_counts);

% Guarantee sparsity
A=sparse(A);
C=column_numbers(:);
R=replication_counts(:);
n=size(A,2);

% Sort columns left to right
[c_sorted,ord]=sort(C);
r_sorted=R(ord);

% Number of extra columns after each original column
extra=zeros(1,n);
extra(c_sorted)=r_sorted-1;

% Shift map for column insertion
shift=[0 cumsum(extra(1:end-1))];

% Original matrix information
m=size(A,1);
[ri,ci,vi]=find(A);
nnz_per_col=accumarray(ci,1,[n 1],@sum,0);
total_new_nnz=nnz(A)+sum((r_sorted-1).*nnz_per_col(c_sorted));

% Preallocate output triplets
ro=zeros(total_new_nnz,1,'like',ri);
co=zeros(total_new_nnz,1,'like',ci);
vo=zeros(total_new_nnz,1,'like',vi);

% Copy non-replicated columns
is_rep=false(n,1); is_rep(c_sorted)=true;
nonrep_idx=~is_rep(ci);

cnt=sum(nonrep_idx);
ro(1:cnt)=ri(nonrep_idx);
co(1:cnt)=ci(nonrep_idx)+shift(ci(nonrep_idx))';
vo(1:cnt)=vi(nonrep_idx);
k=cnt;

% Replicate requested columns
for p=1:numel(c_sorted)
    col=c_sorted(p);
    rep=r_sorted(p);

    mask=ci==col;
    rdup=ri(mask);
    vdup=vi(mask);
    nnzc=numel(rdup);

    base=col+shift(col);

    for j=0:rep-1
        ro(k+1:k+nnzc)=rdup;
        co(k+1:k+nnzc)=base+j;
        vo(k+1:k+nnzc)=vdup;
        k=k+nnzc;
    end
end

% Build output matrix
B=sparse(ro,co,vo,m,n+sum(R-1));

end

% Consistency enforcement
function grumble(A,column_numbers,replication_counts)
if (~isnumeric(A))||(~issparse(A))
    error('A must be a sparse numeric matrix.');
end
if (~isnumeric(column_numbers))||(~isvector(column_numbers))||...
   any(column_numbers<1)||any(mod(column_numbers,1)~=0)
    error('column_numbers must be a vector of positive integers.');
end
if (~isnumeric(replication_counts))||(~isvector(replication_counts))||...
   any(replication_counts<1)||any(mod(replication_counts,1)~=0)||...
   (numel(replication_counts)~=numel(column_numbers))
    error('replication_counts must match column_numbers in size and contain positive integers.');
end
if numel(unique(column_numbers))~=numel(column_numbers)
    error('column_numbers must not contain duplicates.');
end
if any(column_numbers>size(A,2))
    error('column_numbers contain indices outside the matrix.');
end
end

% Капиталисты сами продадут нам верёвку, 
% на которой мы их повесим.
%
% Владимир Ленин

