% Replicates specified columns of a matrix or cell array a 
% specified number of times. Syntax:
%
%            B=repcols(A,col_nums,rep_counts)
%
% Parameters:
%
%   A          - a numeric matrix or a cell array
%
%   col_nums   - vector of column indices to replicate
%
%   rep_counts - vector of positive integers specifying
%                how many copies of each column to make
%
% Output:
%
%   B          - same type as A
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=repcols.m>

function B=repcols(A,col_nums,rep_counts)

% Check consistency
grumble(A,col_nums,rep_counts);

% Replication counts for every column
n=size(A,2); rep_map=ones(1,n); 
rep_map(col_nums)=rep_counts(:);

% Build column index vector
col_idx=repelem(1:n,rep_map);

% Extract and replicate
B=A(:,col_idx);

end

% Consistency enforcement
function grumble(A,col_nums,rep_counts)
if (~isnumeric(A))&&(~iscell(A))
    error('A must be numeric or a cell array.');
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
