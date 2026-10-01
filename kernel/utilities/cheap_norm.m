% The cheapest norm for various representations of matrices. CUDA
% stores matrices by rows, Matlab by columns, and polyadic objects
% can only multiply vectors, in which case Algorithm 2.4 from Hig-
% ham and Tisseur's paper:
%
%            https://doi.org/10.1137/S0895479899356080
%
% is used. Syntax:
%
%                    n=cheap_norm(A,t,itmax)
%
% Parameters:
%
%     A     - a matrix, or a polyadic representation thereof
%
%     t     - (optional) number of probe columns in the poly-
%             adic norm estimator, defaults to 1
%
%     itmax - (optional) maximum number of estimator iterati-
%             ons, defaults to 5
%
% Outputs:
%
%     n - infinity-norm for GPU arrays, 1-norm for CPU arrays,
%         and a lower-bound 1-norm estimate for polyadics
%
% Note: some norms are vastly more expensive than others, this
%       function uses the cheapest ones available.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cheap_norm.m>

function n=cheap_norm(A,t,itmax)

% Set the defaults
if nargin<2, t=1;     end
if nargin<3, itmax=5; end

% Check consistency
grumble(A,t,itmax);

% Inf-norm is best on GPU
if isa(A,'gpuArray')
    n=norm(A,inf); return
end

% 1-norm is best on CPU
if ~isa(A,'polyadic')
    n=norm(A,1); return
end

% Relevant dimensions
row_dim=size(A,1); col_dim=size(A,2);

% Respect small dimensions
t=min(t,col_dim);

% Starting matrix
X=ones(col_dim,t);
if t>1
    for col_idx=2:t
        X(:,col_idx)=2*randi([0 1],col_dim,1)-1;
    end
    for col_idx=1:t
        sign_pool=X(:,[1:(col_idx-1) (col_idx+1):t]);
        ntries=0;
        while any(abs(sign_pool'*X(:,col_idx))==col_dim)
            X(:,col_idx)=2*randi([0 1],col_dim,1)-1;
            sign_pool=X(:,[1:(col_idx-1) (col_idx+1):t]);
            ntries=ntries+1;
            if ntries>100
                break
            end
        end
    end
end
X=X/col_dim;

% Estimator state
idx_hist=[]; idx=1:t; idx_best=1; est_old=0; S=zeros(row_dim,t);

% Iteration loop
for k=1:(itmax+1)

    % Matrix-vector products with the operator
    Y=A*X; col_norms=sum(abs(Y),1);

    % Extract the current estimate
    [est,best_col]=max(col_norms);
    if (est>est_old)||(k==2)
        if k>=2, idx_best=idx(best_col); end
    end

    % Stop if the estimate has not improved
    if (k>=2)&&(est<=est_old)
        n=est_old; return
    end

    % Keep track of the best estimate
    est_old=est; S_old=S;
    if k>itmax
        n=est_old; return
    end

    % Phase matrix with the zero convention from the paper
    S=sign(Y); S(S==0)=1;

    % Remove redundant sign probes in the real case
    if isreal(A)
        if all(any(abs(S_old'*S)==row_dim,1))
            n=est_old; return
        end
        if t>1
            for col_idx=1:t
                sign_pool=[S(:,[1:(col_idx-1) (col_idx+1):t]) S_old];
                ntries=0;
                while any(abs(sign_pool'*S(:,col_idx))==row_dim)
                    S(:,col_idx)=2*randi([0 1],row_dim,1)-1;
                    sign_pool=[S(:,[1:(col_idx-1) (col_idx+1):t]) S_old];
                    ntries=ntries+1;
                    if ntries>100
                        n=est_old; return
                    end
                end
            end
        end
    end

    % Adjoint products with the phase matrix
    Z=A'*S; row_scores=max(abs(Z),[],2);

    % Stop when the best column has been reached
    if (k>=2)&&(max(row_scores)==row_scores(idx_best))
        n=est_old; return
    end

    % Follow the largest row scores
    [~,idx]=sort(row_scores,'descend'); idx=idx(:).';

    % Avoid repeated unit-vector probes
    if t>1
        if all(ismember(idx(1:t),idx_hist))
            n=est_old; return
        end
        idx=[idx(~ismember(idx,idx_hist)) idx(ismember(idx,idx_hist))];
    end

    % Assemble the next probe matrix
    X=zeros(col_dim,t);
    for col_idx=1:t
        X(idx(col_idx),col_idx)=1;
    end
    idx_hist=unique([idx_hist idx(1:t)],'stable');

end

% Return the best lower bound encountered
n=est_old;

end

% Consistency enforcement
function grumble(A,t,itmax)
if ~isnumeric(A)
    error('A must be a matrix or a polyadic.');
end
if (~isnumeric(t))||(~isreal(t))||(~isscalar(t))||(t<1)||...
   (mod(t,1)~=0)
    error('t must be a positive real integer.');
end
if (~isnumeric(itmax))||(~isreal(itmax))||(~isscalar(itmax))||...
   (itmax<2)||(mod(itmax,1)~=0)
    error('itmax must be a real integer greater than one.');
end
end

% "Phobia" implies that these misgivings are irrational,
% when they are anything but.
%
% Rod Liddle

