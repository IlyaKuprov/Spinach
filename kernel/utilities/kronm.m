% Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening 
% Kronecker products. Syntax:
%
%                      x=kronm(Q,x)
%
% Parameters:
%
%    Q -    cell array of Kronecker terms
%
%    x -    a vector or a matrix of appropriate dimension
%
% Output:   
%
%    x -    a vector or a matrix of appropriate dimension
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kronm.m>

function x=kronm(Q,x)

% Check consistency
grumble(Q,x);

% Number of matrices in Q
nmats=numel(Q);

% Number of columns in x
ncols=size(x,2);

% Row and column counts in Q
row_dims=zeros(1,nmats); 
col_dims=zeros(1,nmats);
for n=1:nmats
    [row_dims(n),col_dims(n)]=size(Q{nmats-n+1});
end

% Dimension map for x
x_dims=[col_dims,ncols];

% Reshape into the map
x=reshape(full(x),x_dims);

% Run the products
for n=1:nmats
    
    % Shortcut for opium objects
    if isa(Q{nmats-n+1},'opium')&&(Q{nmats-n+1}.coeff~=1)
        x=Q{nmats-n+1}.coeff*x; continue
    end
    
    % If our dimension is first, do not permute
    if n==1
        
        % Unroll other dimensions
        x=reshape(x,[x_dims(1) prod(x_dims)/x_dims(1)]);
        
        % Run multiplication and update dimension map
        x=Q{nmats}*x; x_dims(1)=row_dims(1);
        
        % Roll other dimensions back up
        x=reshape(full(x),x_dims);
        
    % Otherwise, permute() is unavoidable
    else
        
        % Bring forward n-th dimension
        dims=1:numel(x_dims); dims(n)=[];
        dims=[n,dims]; x=permute(x,dims); %#ok<AGROW>
        
        % Unroll other dimensions
        x=reshape(x,[col_dims(n),numel(x)/col_dims(n)]);
        
        % Run multiplication and update dimension map
        x=Q{nmats-n+1}*x; x_dims(n)=row_dims(n);
        
        % Roll other dimensions back up
        x=reshape(full(x),[row_dims(n),x_dims(dims(2:end))]);
        
        % Put the current dimension back
        x=ipermute(x,dims);
        
    end
    
end

% Return to output dimensions
x=reshape(x,[prod(row_dims),ncols]);

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

