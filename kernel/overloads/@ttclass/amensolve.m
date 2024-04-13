% Solves the linear system Ax=y using the AMEn iteration. Syntax:
%  
%                   x=amensolve(A,y,tol,opts,x0)
%
% Parameters:  
%
%   A    - ttclass representing a square matrix
%
%   y    - ttclass representing the right hand side
%
%   tol  - relative approximation and stopping tolerance,
%          1e-6 is a good start
%
%   x0   - ttclass representing the initial guess
%
% Options (pass empty array for defaults):
%
%   opts.nswp - maximum number of AMEn sweeps
%
%   opts.init_guess_rank - the rank of the initial guess
%
%   opts.enrichment_rank - the rank of the residual and 
%                          enrichment
%
%   opts.resid_damp - local accuracy gap
%
%   opts.rmax - maximum TT rank limit for the solution
%
%   opts.max_full_size - direct vs iterative solver switch-
%                        over dimension
%
%   opts.local_iters - maximum number of bicgstab iterations
%                      for local problems
%
%   opts.verb - Verbosity level: silent (0), sweep (1) or full (2)
%
% Outputs: 
%
%   x    - ttclass representing the solution such that 
%          |x-X| < tol |X| in Frobenius norm, where X is
%          the exact solution.
%
% Note: A, y and x0 should have ntrains==1. Call shrink()
%       on all three if that is not the case.
%
% d.savosyanov@soton.ac.uk
% sergey.v.dolgov@gmail.com
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/amensolve.m>

function x=amensolve(A,y,tol,opts,x0)

% Get the default options started
if isempty(opts), opts=struct; end

% Maximum number of AMEn sweeps
if ~isfield(opts, 'nswp');            opts.nswp=20;           end

% The rank of the initial guess
if ~isfield(opts, 'init_guess_rank'); opts.init_guess_rank=2; end

% The rank of the residual and enrichment
if ~isfield(opts, 'enrichment_rank'); opts.enrichment_rank=4; end

% Local accuracy gap
if ~isfield(opts, 'resid_damp');      opts.resid_damp=2;      end

% Maximum TT rank limit for the solution
if ~isfield(opts, 'rmax');            opts.rmax=Inf;          end

% Direct vs iterative solver switchover dimension
if ~isfield(opts, 'max_full_size');   opts.max_full_size=500; end

% Maximum number of bicgstab iterations for local problems
if ~isfield(opts, 'local_iters');     opts.local_iters=100;   end

% Verbosity level: silent (0), sweep (1) or full (2)
if ~isfield(opts, 'verb');            opts.verb=1;            end

% Initial guess
if ~exist('x0','var')
    [x0,~]=ttort(rand(y,opts.init_guess_rank),-1);
end
x=x0;

% Check inputs for consistency
grumble(A,y,x);

% Read dimension and mode sizes
d=y.ncores;     % dimension of the problem
sz=A.sizes;     % mode sizes
sz=sz(:,1);     % square matrix is assumed

% Clear coefficients from all data
A=clearcoeff(A);
y=clearcoeff(y);
x=clearcoeff(x);

% Initialize reductions of the system to the interfaces of the solution
reduction_XAX = cell(d+1,1); reduction_XAX{1}=1; reduction_XAX{d+1}=1;
reduction_XY = cell(d+1,1);  reduction_XY{1}=1; reduction_XY{d+1}=1;

% If the residual is used, initialise its ranks and reductions
if opts.enrichment_rank>0    
    
    % Residual ranks
    rz = [1;opts.enrichment_rank*ones(d-1,1);1];
    
    % Reductions of the system to the interfaces of the residual
    reduction_ZAX = cell(d+1,1); reduction_ZAX{1}=1; reduction_ZAX{d+1}=1;
    reduction_ZY = cell(d+1,1); reduction_ZY{1}=1; reduction_ZY{d+1}=1;
    
end

% The TT truncation accumulates the error from a single step with the
% factor sqrt(d). Therefore, the tolerance should be adjusted to sqrt(d).
local_tolerance = tol/sqrt(d);

% Initialize some variables of the main loop
flipped = false;          % It will mark that the TT formats are flipped
satisfied = false;        % Check whether error or iteration stopping condition is hit
iter = 1;                 % The number of sweeps

% To prevent overflow, we keep the norms of A,y and x separately in the factored form
norm_x = ones(d-1,1);     % norms of solution blocks
norm_yAx = ones(d-1,1);   % Rescaling factors, equal to nrm(y(i))/nrm(A(i))/nrm(x(i))

% AMEn iteration, loop over sweeps
while ~satisfied        
    
    % Read TT ranks
    ra=A.ranks; ry=y.ranks; rx=x.ranks;
    maximal_error=0;          % Maximal relative error 
    maximal_residual=0;       % Max residual over the sweep
    
    % Loop over cores
    for k=d:-1:1
        
        if ~flipped
            
            % Iterating from d to 1, just orthogonalize the blocks
            current_block = reshape(x.cores{k,1}, [rx(k), sz(k)*rx(k+1)]);
            [v,u]=qr(current_block.',0); new_rx=size(v,2); u=u.'; v=v.';
            
        else
            
            % Assemble the local right-hand side
            current_rhs = local_vector(reduction_XY{k}, y.cores{k,1}, reduction_XY{k+1});
            current_rhs = current_rhs*prod(norm_yAx);
            
            % Read the initial guess -- the TT block from the prev. iteration
            old_block = reshape(x.cores{k,1}, [rx(k)*sz(k)*rx(k+1),1]);
            
            % Save the norm of the right-hand side
            norm_rhs = norm(current_rhs,2);
            
            % Read the matrix parts, accelerate a plenty of iterations with them
            left_XAX = reduction_XAX{k};
            current_A = A.cores{k,1};
            right_XAX = reduction_XAX{k+1};
            
            % Measure the previous residual
            previous_Ax=local_matvec(old_block, left_XAX, current_A, right_XAX);
            previous_residual = norm(current_rhs-previous_Ax,2)/norm_rhs;
            
            % Look at the size
            if (rx(k)*sz(k)*rx(k+1)<opts.max_full_size)
                
                % Assemble the full system and solve directly
                current_matrix = local_matrix(left_XAX, current_A, right_XAX);
                current_block = current_matrix\current_rhs;
                
                % Report if verbocity is >1
                if (opts.verb>1)
                    fprintf('amen_solve: swp=%d, i=%d. Backslash. ', iter, k);
                end
                
            else
                
                % Extract the number of local iterations wisely                
                local_iters = opts.local_iters;
                if (numel(local_iters)>1)
                    local_iters = local_iters(k);
                end
                
                % Run the bicgstab with the matrix defined by the structured MatVec
                [current_block,~,bicg_res,bicg_iters] = bicgstab(...
                    @(x)local_matvec(x, left_XAX, current_A, right_XAX),...
                    current_rhs, max(local_tolerance/opts.resid_damp,eps*2),...
                    local_iters, [], [], old_block);
                
                % Report to user
                if (opts.verb>1)
                    fprintf('amen_solve: iter=%d, i=%d. Bicgstab: %g iters, residual %3.3e. ', iter, k, bicg_iters, bicg_res);
                end
                
            end
            
            % Measure the error
            local_error = norm(current_block-old_block,2)/norm(current_block,2);
            
            % Update the worst errors
            maximal_error = max(maximal_error, local_error);
            maximal_residual = max(maximal_residual, previous_residual);
            
            % Reshape the block for the truncation
            current_block = reshape(current_block, [rx(k), sz(k)*rx(k+1)]);  
            
            % Truncation
            if k>1
                
                % Compute the SVD
                [u,s,v] = svd(current_block,'econ'); s = diag(s);
                
                % Select the rank based on Fro-norm thresholding
                new_rx = frob_chop(s,local_tolerance*norm(s,2)); 
                
                % Limit the rank to rmax
                new_rx = min(new_rx, opts.rmax);
                
                % Shrink SVD factors to r
                % U*S will be cast to the neighbouring block
                u = u(:,1:new_rx)*diag(s(1:new_rx));
                
                % The current block will be left orthogonal
                v = v(:,1:new_rx)';
                
                % Save the truncated solution to compute the residual later
                current_block = u*v;
                
                % Report the chosen rank
                if (opts.verb>1)
                    fprintf('Rank: %d. \n', new_rx);
                end
                
            end
            
            % Copy the new block to the main storage
            x.cores{k} = reshape(current_block, rx(k), sz(k), 1, rx(k+1));
            
        end
        
        % Update the residual z
        if (opts.enrichment_rank>0)
            
            % In the first sweep, there are no reductions yet
            if iter==1
                
                % Just initialize z by random
                zblock = randn(rz(k), sz(k)*rz(k+1));
                
            % In higher sweeps, the reductions are computed
            else
                
                % Project y via the reductions onto the Z interface
                zblock_y = local_vector(reduction_ZY{k}, y.cores{k,1}, reduction_ZY{k+1});
                zblock_y = zblock_y*prod(norm_yAx);
                
                % Project Ax
                zblock_Ax = local_matvec(x.cores{k,1}, reduction_ZAX{k}, A.cores{k,1}, reduction_ZAX{k+1});
                
                % z=y-Ax projected to Z
                zblock = zblock_y-zblock_Ax;
                
            end
            
            % Reshape zblock for the orthogonalization
            zblock = reshape(zblock, [rz(k), sz(k)*rz(k+1)]);
            [zblock,~]=qr(zblock.', 0); zblock = zblock.';
            
            % Careful: store the old rank of z, since it is that will be used
            % in the solution enrichment, not the updated value after the QR
            old_rz = rz(k);
            
            % Now replace it
            rz(k) = size(zblock,1);
            zblock = reshape(zblock, [rz(k), sz(k), 1, rz(k+1)]);
            
        end
        
        % Enrichment
        if k>1
            
            % Apply enrichment to the solution
            if (opts.enrichment_rank>0)&&(iter>1)&&(flipped)
                
                % The enrichment uses a mixed interface: the one of Z on
                % the left, but the one of X on the right, since it must
                % concatenate with the solution block.
                
                % First, project y to Z-I-X
                zblock_y = local_vector(reduction_ZY{k}, y.cores{k,1}, reduction_XY{k+1});
                zblock_y = zblock_y*prod(norm_yAx);
                
                % Project Ax to Z-I-X
                zblock_Ax = local_matvec(x.cores{k,1}, reduction_ZAX{k}, A.cores{k,1}, reduction_XAX{k+1});
                
                % z=y-Ax projected to Z-I-X
                enrichment_block = zblock_y-zblock_Ax;
                enrichment_block = reshape(enrichment_block, [old_rz, sz(k)*rx(k+1)]);

                % Enrichment is made here
                v2 = [v; enrichment_block]; v2 = v2.';
                [v,rv]=qr(v2, 0); v = v.';
                
                % Pass the L-factor to the next block
                rv = rv(:,1:new_rx); u = u*rv.';
                
                % Update the current solution rank
                new_rx = size(v,1);
                
            end
            
            % Measure and remove the norm of u
            new_norm_x = norm(u,'fro');
            if (new_norm_x>0)
                u = u/new_norm_x;
            else
                new_norm_x=1;
            end
            
            % Store it in norm_x
            norm_x(k-1)=norm_x(k-1)*new_norm_x;
            
            % Calculate the White's prediction for the next block by casting u onto it
            next_block = x.cores{k-1,1};
            next_block = reshape(next_block, [rx(k-1)*sz(k-1), rx(k)]);
            next_block = next_block*u;
            
            % Copy the rank
            rx(k) = new_rx;
            
            % Now it is a good initial guess
            x.cores{k-1,1} = reshape(next_block, rx(k-1), sz(k-1),1, rx(k)); 
            
            % Copy the orthogonal current block to the storage
            v = reshape(v, [rx(k), sz(k),1, rx(k+1)]);
            x.cores{k,1} = v;
            
            % Compute next interface reductions with X
            [reduction_XAX{k},norm_A] = reduce_matrix(reduction_XAX{k+1}, v, A.cores{k,1}, v);
            [reduction_XY{k},norm_y] = reduce_vector(reduction_XY{k+1}, v, y.cores{k,1});
            
            % Compute next interface reductions with Z
            if (opts.enrichment_rank>0)
                
                % Run the reductions
                reduction_ZAX{k} = reduce_matrix(reduction_ZAX{k+1}, zblock, A.cores{k,1}, v);
                reduction_ZY{k} = reduce_vector(reduction_ZY{k+1}, zblock, y.cores{k,1});
                
                % Remove the norm of A and y
                reduction_ZAX{k} = reduction_ZAX{k}/norm_A;
                reduction_ZY{k} = reduction_ZY{k}/norm_y;
                
            end
            
            % Compute the new rescaling factor
            norm_yAx(k-1) = (norm_y/(norm_A*norm_x(k-1)));
            
        end
        
    end
    
    % Flip reductions and the ranks of Z
    reduction_XAX = revert_interface(reduction_XAX, ra, rx, rx);
    reduction_XY  = revert_interface(reduction_XY, ry, rx);
    if (opts.enrichment_rank>0)
        reduction_ZAX = revert_interface(reduction_ZAX, ra, rz, rx);
        reduction_ZY  = revert_interface(reduction_ZY, ry, rz);
        rz = rz(d+1:-1:1, :);
    end
    
    % Flip tensor trains
    y=revert(y); A=revert(A); x=revert(x);
    
    % Flip sizes and norms
    sz=sz(d:-1:1); norm_x = norm_x(d-1:-1:1);
    norm_yAx = norm_yAx(d-1:-1:1); flipped = ~flipped;
    
    % Report sweep info
    if (opts.verb>0)&&(~flipped)
        fprintf('amen_solve: iter=%d, err=%3.3e, res=%3.3e, rank=%d\n', iter, maximal_error, maximal_residual, max(rx));
    end
    
    % Check the stopping criteria
    satisfied = ((maximal_error<tol)||(iter>=opts.nswp))&&(~flipped);
    
    % Count the number of sweeps
    if ~flipped, iter = iter+1; end
    
end

% Cast spatial solution to the desired form
rx = x.ranks;
norm_x = exp(sum(log(norm_x))/d); % distribute norms equally
for k=1:d
    x.cores{k,1} = reshape(x.cores{k,1}, [rx(k), sz(k), 1, rx(k+1)]); % store the sizes in
    x.cores{k,1} = x.cores{k,1}*norm_x;
end
x.coeff=1; x.tolerance=(norm_x^d)*tol;

end

% Consistency enforcement
function grumble(A,y,x0)
if ~isa(A,'ttclass') || ~isa(y,'ttclass')
    error('Both matrix and right-hand-side should be ttclass.');
end
if A.ncores~=y.ncores
    error('Dimensions of matrix and right-hand-side do not match.');
end
if A.ntrains>1 || y.ntrains>1
    error('Both matrix and right-hand-side should be shrinked before solution.');
end
d=A.ntrains; szA=A.sizes;
if ~all(szA(:,1)==szA(:,2))
    error('Matrix should be square.');
end
szy=y.sizes;
if ~all(szy(:,1)==szA(:,2))
    error('Mode sizes of the right-hand-side do not match those of the matrix.')
end
if ~all(szy(:,2)==ones(d,1))
    error('The right-hand-side should be a vector.')
end
if exist('x0','var')
    if A.ncores~=x0.ncores
        error('Dimensions of matrix and initial guess do not match.');
    end
    szx=x0.sizes;
    if ~all(szx(:,1)==szA(:,2))
        error('Mode sizes of the initial guess do not match those of the matrix.')
    end
    if ~all(szx(:,2)==ones(d,1))
        error('The initial guess should be a vector.')
    end
end
end

% Takes right reductions of sizes (ra x rw x rx), and flips them to the left
% form (rw x rx x ra). Note: ranks should be given in the initial order.
function [interface] = revert_interface(interface, ra, rw, rx)
d = numel(interface)-1; % We trust it to come with R==1
if (nargin<4)
    rx = ones(size(rw));
end
for i=1:d
    interface{i} = reshape(interface{i}, ra(i), rw(i)*rx(i));
    interface{i} = interface{i}.';
    if (nargin>3)
        interface{i} = reshape(interface{i}, rw(i), rx(i), ra(i));
    end
end
interface = interface(d+1:-1:1);
end

% Accumulates the right interface reduction of the matrix, W{k:d}'*A{k:d}*X{k:d}
% Right reduction has the form of the last matrix TT block, i.e. [ra, rw, rx]
function [interface,nrm] = reduce_matrix(interface, w, A, x)
[ra1,n,m,ra2]=size(A);
[rx1,~,~,rx2]=size(x);
[rw1,~,~,rw2]=size(w);
wc = reshape(w, rw1, n*rw2);
wc = conj(wc);
xc = reshape(x, rx1*m, rx2);
interface = reshape(interface, ra2*rw2, rx2);
interface = xc*interface.'; % size rx1 m x ra2 rw2
interface = reshape(interface, rx1, m*ra2*rw2);
interface = interface.';
interface = reshape(interface, m*ra2, rw2*rx1);
tmp = reshape(A, ra1*n, m*ra2);
interface = tmp*interface;  % size ra1(k)*n, rw2*rx1
interface = reshape(interface, ra1, n*rw2*rx1);
interface = interface.';
interface = reshape(interface, n*rw2, rx1*ra1);
interface = wc*interface;   % size rw1, rx1 ra1
interface = reshape(interface, rw1*rx1, ra1);
interface = interface.';
if (nargout>1)
    nrm = norm(interface,'fro');
    interface = interface/nrm;
end
interface = reshape(interface, ra1, rw1, rx1);
end

% Accumulates the right interface reduction of the vector, W{k:d}'*X{k:d}
% Right reduction has the form of the last vector TT block, i.e. [rx, rw]
function [interface,nrm] = reduce_vector(interface, w, x)
[rw1,n,~,rw2]=size(w);
[rx1,~,~,rx2]=size(x);
wc = reshape(w, rw1, n*rw2);
tmp = reshape(x, rx1*n, rx2);
interface = tmp*interface; % size rx1 n x rw2
interface = reshape(interface, rx1, n*rw2);
interface = interface*wc'; % size rx1, rw1
if (nargout>1)
    nrm = norm(interface,'fro');
    interface = interface/nrm;
end
end

% A matrix-vectors product for the matrix in the 3D TT (WAX1-A-WAX2), and
% full vectors of size (rx1*m*rx2). Returns (rw1*n*rw2)
function w=local_matvec(x, left_WAX, A, right_WAX)
[ra1,n,m,ra2]=size(A);
[rw1,rx1,~] = size(left_WAX);
[~,rw2,rx2] = size(right_WAX);
xc = reshape(x, [rx1*m, rx2]);
tmp = reshape(right_WAX, [ra2*rw2, rx2]);
w = xc*tmp.';
w = reshape(w, rx1, m*ra2*rw2);
w = w.';
w = reshape(w, m*ra2, rw2*rx1);
tmp = reshape(A, ra1*n, m*ra2);
w = tmp*w;
w = reshape(w, ra1*n*rw2, rx1);
w = w.';
w = reshape(w, rx1*ra1, n*rw2);
tmp = reshape(left_WAX, rw1, rx1*ra1);
w = tmp*w;
w = reshape(w, rw1*n*rw2, 1);
end

% Builds the full (rw1*n*rw2) x (rx1*m*rx2) matrix from its TT blocks
function B=local_matrix(left_WAX, A, right_WAX)
[ra1,n,m,ra2]=size(A);
[rw1,rx1,~] = size(left_WAX);
[~,rw2,rx2] = size(right_WAX);
B = reshape(left_WAX, [rw1*rx1, ra1]);
tmp = reshape(A, [ra1, n*m*ra2]);
B = B*tmp;
B = reshape(B, [rw1, rx1, n, m*ra2]);
B = permute(B, [1,3,2,4]);
B = reshape(B, [rw1*n*rx1*m, ra2]);
tmp = reshape(right_WAX, [ra2, rw2*rx2]);
B = B*tmp;
B = reshape(B, [rw1*n, rx1*m, rw2, rx2]);
B = permute(B, [1,3,2,4]);
B = reshape(B, [rw1*n*rw2, rx1*m*rx2]);
end

% Builds the full (rw1*n*rw2) x 1 vector from its TT blocks
function w=local_vector(left_WX, x, right_WX)
[rx1,n,~,rx2]=size(x);
[rw1,~] = size(left_WX);
[~,rw2] = size(right_WX);
w = reshape(x, [rx1, n*rx2]);
w = left_WX*w;
w = reshape(w, [rw1*n, rx2]);
w = w*right_WX;
w = reshape(w, [rw1*n*rw2, 1]);
end

% Had I the heavens' embroidered cloths,
% Enwrought with golden and silver light,
% The blue and the dim and the dark cloths
% Of night and light and the half light,
% I would spread the cloths under your feet:
% But I, being poor, have only my dreams;
% I have spread my dreams under your feet;
% Tread softly because you tread on my dreams.
%
% W.B. Yeats (1865-1939)

