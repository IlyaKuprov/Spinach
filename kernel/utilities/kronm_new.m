% Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening 
% Kronecker products. Syntax:
%
%                      M=kronm(Q,M)
%
% Parameters:
%
%    Q -    cell array of Kronecker terms
%
%    M -    a vector or a matrix of appropriate dimension
%
% Output:   
%
%    M -    a vector or a matrix of appropriate dimension
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kronm.m>

function M=kronm_new(Q,M)

% Check consistency
grumble(Q,M);

% Dimension statistics
n_mats_in_q=numel(Q);
n_cols_in_m=size(M,2);

% Row and column counts in Q
row_dims=zeros(1,n_mats_in_q); 
col_dims=zeros(1,n_mats_in_q);
for n=1:n_mats_in_q
    [row_dims(n),col_dims(n)]=size(Q{n_mats_in_q-n+1});
end

% Fold up implicit dimensions of M
M=reshape(full(M),[col_dims n_cols_in_m]);

% Run the products
for n=1:n_mats_in_q

    % Contract each implicit dimension
    M=tensorprod(full(Q{n}),M,2,n_mats_in_q);
    
end

% Unfold implicit dimensions of M
M=reshape(M,[prod(row_dims) n_cols_in_m]);

end

% Consistency enforcement
function grumble(Q,x)
if (~iscell(Q))
    error('Q must be a cell array.');
end
for n=1:numel(Q)
    if ~ismatrix(Q{n})
        error('Q must be a cell array of matrices.');
    end
end
if ~isnumeric(x)
    error('x must be numeric.');
end
end

% Never get upset at anyone. Either forgive
% a man, or kill him.
%
% Joseph Stalin

