% The cheapest possible norm for various representations of matrices. CUDA 
% stores matrices by rows, Matlab by columns, and polyadic objects can on-
% ly multiply vectors, in which case Algorithm 2.1 from Hager's paper:
%
%                 https://doi.org/10.1137/S0895479899356080
%
% is used with Higham's extension to the complex case. Syntax:
%
%                             n=cheap_norm(A)
%
% Parameters:
%
%     A - a matrix, or a polyadic representation thereof
%
% Outputs:
%
%     n - infinity-norm for GPU arrays, 1-norm for CPU arrays,
%         and Hager's upper bound on 1-norm for polyadics
%
% Note: some norms are vastly more expensive than others, this function
%       uses the cheapest ones available.
%
% Note: in very rare cases, Hager's algorithm stagnates. If this happens,
%       drop us a note and we will implement Higham's more recent work.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cheap_norm.m>

function n=cheap_norm(A)

% Check consistency
grumble(A);

% Inf-norm is best on GPU
if isa(A,'gpuArray')
    n=norm(A,inf); return
end

% 1-norm is best on CPU
if ~isa(A,'polyadic')
    n=norm(A,1); return
end

% Relevant dimension
dim=size(A,2);

% Starting vector
x=ones(dim,1)/dim;

% Column index and counter
idx_prev=0; n_iter=0;

% Iteration loop
while true
    
    % Hager's magic
    y=A*x; z=A'*sign(y);
    [~,idx]=max(abs(z));
    
    % Termination condition
    if idx==idx_prev
        n=norm(y,1); break
    end
    
    % Follow the largest element
    x=zeros(dim,1); x(idx)=1;
    
    % Column index update
    idx_prev=idx;
    
    % Iteration counter
    n_iter=n_iter+1;
    
    % Fallback scenario
    if n_iter>10
        error('norm estimator stagnated.');
    end
    
end

end

% Consistency enforcement
function grumble(A)
if ~isnumeric(A)
    error('A must be a matrix or a polyadic.');
end
end

% "Phobia" implies that these misgivings are irrational,
% when they are anything but.
%
% Rod Liddle

 