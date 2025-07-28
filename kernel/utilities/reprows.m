% Replicates specified rows of a matrix a specified num-
% ber of times. Syntax:
%
%            B=reprows(A,row_nums,rep_counts)
%
% Parameters:
%
%   A          - a sparse matrix
%
%   row_nums   - vector of row indices to replicate
%
%   rep_counts - vector of positive integers specifying
%                how many copies of each row to make
%
% Output:
%
%   B          - a sparse matrix
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=reprows.m>

function B=reprows(A,row_nums,rep_counts)

% Check consistency
grumble(A,row_nums,rep_counts);

% Sort user indices top to bottom
R=row_nums(:); S=rep_counts(:);
[R,ord]=sort(R); S=S(ord);

% Number of extra rows after each original row
n=size(A,1); extra_rows=zeros(1,n); extra_rows(R)=S-1;

% Shift map for row insertion
shift=[0 cumsum(extra_rows(1:end-1))];

% Original matrix information
m=size(A,2); [ri,ci,vi]=find(A);
nnz_per_row=accumarray(ri,1,[n 1],@sum,0);
total_new_nnz=nnz(A)+sum((S-1).*nnz_per_row(R));

% Preallocate output triplets
ro=zeros(total_new_nnz,1,'like',ri);
co=zeros(total_new_nnz,1,'like',ci);
vo=zeros(total_new_nnz,1,'like',vi);

% Copy non-replicated rows
is_rep=false(n,1); is_rep(R)=true;
nonrep_idx=~is_rep(ri);
cnt=sum(nonrep_idx);
ro(1:cnt)=ri(nonrep_idx)+shift(ri(nonrep_idx))';
co(1:cnt)=ci(nonrep_idx);
vo(1:cnt)=vi(nonrep_idx); k=cnt;

% Replicate requested rows
for p=1:numel(R)

    row=R(p); rep=S(p);
    mask=(ri==row);
    cdup=ci(mask);
    vdup=vi(mask);
    nnzr=numel(cdup);
    base=row+shift(row);

    for j=0:(rep-1)
        ro(k+1:k+nnzr)=base+j;
        co(k+1:k+nnzr)=cdup;
        vo(k+1:k+nnzr)=vdup;
        k=k+nnzr;
    end

end

% Build the output matrix
B=sparse(ro,co,vo,n+sum(S-1),m);
if ~issparse(A), B=full(B); end

end

% Consistency enforcement
function grumble(A,row_nums,rep_counts)
if ~isnumeric(A)
    error('A must be a numeric matrix.');
end
if (~isnumeric(row_nums))||(~isvector(row_nums))||...
   any(row_nums<1)||any(mod(row_nums,1)~=0)
    error('row_nums must be a vector of positive integers.');
end
if (~isnumeric(rep_counts))||(~isvector(rep_counts))||...
   any(rep_counts<1)||any(mod(rep_counts,1)~=0)||...
   (numel(rep_counts)~=numel(row_nums))
    error('rep_counts must match row_nums in size and contain positive integers.');
end
if numel(unique(row_nums))~=numel(row_nums)
    error('row_nums must not contain duplicates.');
end
if any(row_nums>size(A,1))
    error('row_nums indices exceed matrix dimension.');
end
end

% It's important to draw a distinction between
% loneliness and comfortable solitude; the for-
% mer is a tragedy, the latter a blessing.
%
% Hannah Tomes, in The Spectator

