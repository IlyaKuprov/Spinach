% Replicates specified rows of a matrix or cell array a specified
% number of times. Syntax:
%
%            B=reprows(A,row_nums,rep_counts)
%
% Parameters:
%
%   A          - a numeric matrix or a cell array
%
%   row_nums   - vector of row indices to replicate
%
%   rep_counts - vector of positive integers specifying
%                how many copies of each row to make
%
% Output:
%
%   B          - same type as A
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=reprows.m>

function B=reprows(A,row_nums,rep_counts)

% Check consistency
grumble(A,row_nums,rep_counts);

% Replication counts for every row
n=size(A,1); rep_map=ones(1,n); rep_map(row_nums)=rep_counts(:);

% Build row index vector
row_idx=repelem(1:n,rep_map);

% Extract and replicate
B=A(row_idx,:);

end

% Consistency enforcement
function grumble(A,row_nums,rep_counts)
if (~isnumeric(A))&&(~iscell(A))
    error('A must be numeric or a cell array.');
end
if ndims(A)~=2
    error('A must be two-dimensional.');
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
