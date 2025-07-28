% REPCOLS  Replicate selected sparse columns in‑place.
%
% B = REPCOLS(A, C, R) returns a sparse matrix B obtained from the sparse
% matrix A by replacing each column C(i) with R(i) identical consecutive
% copies.  All other columns stay in the same left‑to‑right order.
%
% INPUT:
%   A  – m‑by‑n sparse matrix (any numeric class)
%   C  – vector of 1‑based column indices to be replicated
%   R  – vector of positive integers, same size as C, telling how many
%        copies to make of the corresponding column in C
%
% OUTPUT:
%   B  – m‑by‑(n + sum(R‑1)) sparse matrix
%
% EXAMPLE:
%   A = sparse([1 0 2 ; 0 3 0 ; 4 0 5]);
%   % replicate column 1 three times, column 3 twice:
%   B = repcols(A, [1 3], [3 2]);
%
% ChatGPT o3
% ilya.kuprov@weizmann.ac.il

function B = repcols(A, column_numbers, replication_counts)

    %--------------------------- validation ------------------------------%
    if nargin~=3
        error('Exactly three inputs required: B = repcols(A,C,R).');
    end
    A  = sparse(A);                     % guarantees sparsity
    C  = column_numbers(:);
    R  = replication_counts(:);
    if ~isequal(size(C),size(R))
        error('column_numbers and replication_counts must be the same size.');
    end
    if any(R < 1 | R ~= round(R))
        error('replication_counts must be positive integers.');
    end
    n = size(A,2);
    if any(C < 1 | C > n)
        error('column_numbers contain indices outside 1..size(A,2).');
    end
    if numel(unique(C)) ~= numel(C)
        error('column_numbers must not contain duplicates.');
    end

    %---------------------- preparation & bookkeeping --------------------%
    [Csort,ord] = sort(C);              % process left‑to‑right for clarity
    Rsort       = R(ord);

    extra          = zeros(1,n);        % how many *new* columns after each old
    extra(Csort)   = Rsort - 1;         % e.g. replicate 4 → add 3 extras

    % shift(k) = how many columns have already been *inserted* before col k
    shift          = [0 cumsum(extra(1:end-1))];

    m  = size(A,1);
    [ri,ci,vi]     = find(A);           % original triplets
    nnz_per_col    = accumarray(ci,1,[n 1],@sum,0);
    total_new_nnz  = nnz(A) + sum((Rsort-1).*nnz_per_col(Csort));

    % preallocate output triplets
    ro = zeros(total_new_nnz,1,'like',ri);
    co = zeros(total_new_nnz,1,'like',ci);
    vo = zeros(total_new_nnz,1,'like',vi);

    %---------------------- copy non‑replicated columns -------------------%
    is_rep     = false(n,1);  is_rep(Csort) = true;
    nonrep_idx = ~is_rep(ci);           % logical mask on entries, not cols

    cnt        = sum(nonrep_idx);
    ro(1:cnt)  = ri(nonrep_idx);
    co(1:cnt)  = ci(nonrep_idx) + shift(ci(nonrep_idx))';
    vo(1:cnt)  = vi(nonrep_idx);
    k          = cnt;                   % running index

    %---------------------- replicate requested columns ------------------%
    for p = 1:numel(Csort)
        col = Csort(p);
        rep = Rsort(p);

        mask  = (ci == col);            % entries belonging to this column
        rdup  = ri(mask);
        vdup  = vi(mask);
        nnzc  = numel(rdup);

        base  = col + shift(col);       % where first copy will land

        for j = 0:rep-1
            ro(k+1 : k+nnzc) = rdup;
            co(k+1 : k+nnzc) = base + j;
            vo(k+1 : k+nnzc) = vdup;
            k = k + nnzc;
        end
    end

    %--------------------------- build output ----------------------------%
    B = sparse(ro,co,vo, m, n + sum(R-1));

end

% Капиталисты сами продадут нам верёвку, 
% на которой мы их повесим.
%
% Владимир Ленин

