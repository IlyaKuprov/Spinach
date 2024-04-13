% A sum across all dimensions of a cell array. Syntax:
%
%                       S=totsum(A)
%
% Parameters:
%
%      A - a cell array of numerical objects
%
% Outputs:
%
%      S - the sum of all elements in A
%
% Notes: if all elements of A are sparse, a sparse result
%        will be returned.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cell/totsum.m>

function S=totsum(A)

% Check consistency
grumble(A);

% Check array type
r_u_sparse=cellfun(@issparse,A);
sparse_path=any(r_u_sparse(:));

% Run the addition
if sparse_path
    
    % Run sparse matrix addition
    rows=cell(numel(A),1); 
    cols=cell(numel(A),1);
    vals=cell(numel(A),1);
    for n=1:numel(A)
        [rows{n},cols{n},vals{n}]=find(A{n});
    end
    rows=cell2mat(rows); 
    cols=cell2mat(cols); 
    vals=cell2mat(vals);
    S=sparse(rows,cols,vals,size(A{1},1),size(A{1},2));
    
else
    
    % Run full matrix addition
    S=zeros(size(A{1}));
    for n=1:numel(A)
        S=S+A{n};
    end
    
end

end

% Consistency enforcement
function grumble(A)
r_u_numeric=cellfun(@isnumeric,A);
if ~all(r_u_numeric(:))
    error('all elements of A must be numeric.');
end
end

% The idea that global warming is the most important 
% problem facing the world is total nonsense and is
% doing a lot of harm. 
%
% Freeman Dyson

