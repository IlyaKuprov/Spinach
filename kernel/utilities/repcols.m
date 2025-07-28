% Replicates specified columns of a matrix a specified
% number of times. Syntax:
%
%            B=repcols(A,col_num,rep_counts)
%
% Parameters:
%
%   A          - a sparse matrix
%
%   col_nums   - vector of column indices to replicate
%
%   rep_counts - vector of positive integers specifying
%                how many copies of each column to make
%
% Output:
%
%   B          - a sparse matrix
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=repcols.m>

function B=repcols(A,col_nums,rep_counts)

% Check consistency
grumble(A,col_nums,rep_counts);

% Sort user indices left to right
C=col_nums(:); R=rep_counts(:);
[C,ord]=sort(C); R=R(ord);

% Number of extra columns after each original column
n=size(A,2); extra_cols=zeros(1,n); extra_cols(C)=R-1;

% Shift map for column insertion
shift=[0 cumsum(extra_cols(1:end-1))];

% Original matrix information
m=size(A,1); [ri,ci,vi]=find(A);
nnz_per_col=accumarray(ci,1,[n 1],@sum,0);
total_new_nnz=nnz(A)+sum((R-1).*nnz_per_col(C));

% Preallocate output triplets
ro=zeros(total_new_nnz,1,'like',ri);
co=zeros(total_new_nnz,1,'like',ci);
vo=zeros(total_new_nnz,1,'like',vi);

% Copy non-replicated columns
is_rep=false(n,1); is_rep(C)=true;
nonrep_idx=~is_rep(ci);
cnt=sum(nonrep_idx);
ro(1:cnt)=ri(nonrep_idx);
co(1:cnt)=ci(nonrep_idx)+shift(ci(nonrep_idx))';
vo(1:cnt)=vi(nonrep_idx); k=cnt;

% Replicate requested columns
for p=1:numel(C)
    
    col=C(p); rep=R(p);
    mask=(ci==col);
    rdup=ri(mask);
    vdup=vi(mask);
    nnzc=numel(rdup);
    base=col+shift(col);

    for j=0:(rep-1)
        ro(k+1:k+nnzc)=rdup;
        co(k+1:k+nnzc)=base+j;
        vo(k+1:k+nnzc)=vdup;
        k=k+nnzc;
    end

end

% Build the output matrix
B=sparse(ro,co,vo,m,n+sum(R-1));
if ~issparse(A), B=full(B); end

end

% Consistency enforcement
function grumble(A,col_nums,rep_counts)
if ~isnumeric(A)
    error('A must be a numeric matrix.');
end
if (~isnumeric(col_nums))||(~isvector(col_nums))||...
   any(col_nums<1)||any(mod(col_nums,1)~=0)
    error('col_nums must be a vector of positive integers.');
end
if (~isnumeric(rep_counts))||(~isvector(rep_counts))||...
   any(rep_counts<1)||any(mod(rep_counts,1)~=0)||...
   (numel(rep_counts)~=numel(col_nums))
    error('rep_counts must match col_nums in size and contain positive integers.');
end
if numel(unique(col_nums))~=numel(col_nums)
    error('col_nums must not contain duplicates.');
end
if any(col_nums>size(A,2))
    error('col_nums indices exceed matrix dimension.');
end
end

% Капиталисты сами продадут нам верёвку, 
% на которой мы их повесим.
%
% Владимир Ленин

