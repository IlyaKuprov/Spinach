% Propagation step function. Computes the action by a matrix exponential 
% without compuing that exponential. Supports one-, two-, and three-point
% product quadratures. Syntax:
%
%                   rho=step(spin_system,L,rho,time_step)
%
% Parameters:
%
%      L          -  Liouvillian or Hamiltonian to be used for 
%                    propagation; centre point piecewise-constant
%                    rule if one matrix is supplied, piecewise-
%                    linear rule if two matrices {left, right} 
%                    are supplied, piecewise-quadratic if three
%                    matrices {left, midpoint, right} are given.
%
%                    State-dependent evolution generators are 
%                    supported: if L{1} is a function handle (see
%                    iserstep.m documentation), L{2} is current
%                    time, and L{3} is the method (see iserstep.m
%                    documentation), the problem is routed to a
%                    an appropriate Lie group solver. 
%
%      rho        -  state vector or density matrix
%
%      time_step  -  length of the time step to take
%
% Outputs:
%
%      rho        -  state vector or density matrix
%
% Note: we initially had a faithful implementation of the Krylov process
%       here - subspace, orthogonalisation, projection, etc., but in all
%       our testing it was much inferior to the reordered Taylor process
%       that is currently implemented below.
%
% Note: the peculiar sequence of algebraic operations in the code below
%       is designed to minimise the memory footprint in large cases.
%
% Note: set sys.expmv_backend='auto' to enable heuristic selection of the
%       exponential-action backend in single-vector numeric propagation.
%       The default value is 'default', which uses the native reordered
%       Taylor procedure.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
% a.acharya@soton.ac.uk
% c.musselwhite@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=step.m>

function rho=step(spin_system,L,rho,time_step)

% If L{1} is a function handle, route to iserstep
if iscell(L)&&(numel(L)==3)&&isa(L{1},'function_handle')

    % Call state-dependent evolution generator solver
    rho=iserstep(spin_system,L,rho,time_step); return;

end

% Check consistency
grumble(L,rho,time_step);

% Zero time step shortcut
if time_step==0, return; end

% Do we want to run on GPU?
want_gpu=ismember('gpu',spin_system.sys.enable);

% Is the generator already on GPU?
if isnumeric(L)
    if want_gpu&&(~isa(L,'gpuArray'))
        L=gpuArray(L); 
    end
else
    for n=1:numel(L)
        if want_gpu&&(~isa(L{n},'gpuArray'))
            L{n}=gpuArray(L{n});
        end
    end
end

% Is the state already on GPU?
if isnumeric(rho)
    if want_gpu&&(~isa(rho,'gpuArray'))
        rho=gpuArray(rho); 
    end
else
    for n=1:numel(rho)
        if want_gpu&&(~isa(rho{n},'gpuArray'))
            rho{n}=gpuArray(rho{n});
        end
    end
end

% Is it expm(A)*x (wavefunctions or state vectors)?
expm_times_vec=ismember(spin_system.bas.formalism,{'sphten-liouv',...
                                                   'zeeman-liouv',...
                                                   'zeeman-wavef'});
% expm(A)*rho*expm(-A)
if ~expm_times_vec
    
    % Process Lie generators
    if iscell(L)&&(numel(L)==2)

        % Two-point Lie quadrature
        L=isergen(L{1},[],L{2},time_step);

    elseif iscell(L)&&(numel(L)==3)

        % Three-point Lie quadrature
        L=isergen(L{1},L{2},L{3},time_step);

    end

    % Fast bypass for small density matrices
    if size(L,1)<spin_system.tols.small_matrix

        % Use Matlab's expm
        P=expm(-1i*L*time_step);

        % A cell array
        if iscell(rho)

            % Parallel over cells
            parfor n=1:numel(rho)
                rho{n}=P*rho{n}*P';
            end

        else

            % A martrix
            rho=P*rho*P';

        end

        % Shortcut
        return;

    else

        % Subdivide the time step to ensure monotonic convergence
        norm_gen=cheap_norm(L)*abs(time_step); nsteps=ceil(norm_gen/2);

        % Step checks
        if nsteps>1e4

            % Catch the common mistake of supplying wildly unreasonable |L*dt|
            error('either dt is too long, or L is too big: |L*dt|>1e4, check both.');

        elseif nsteps>100

            % Warn user if too many substeps are needed
            report(spin_system,['WARNING: ' num2str(nsteps)...
                                ' substeps required, consider using evolution() here.']);
        end

        % Check if we have a matrix or
        % a cell array of matrices
        if isnumeric(rho)

            % Encapsulate into a cell array
            return_numeric=true(); rho={rho};

        elseif iscell(rho)

            % Keep the cell array
            return_numeric=false();

        end

        % Decide if parallelisation would be beneficial
        parfor_makes_sense=(size(rho{1},1)>256)&&(numel(rho)>32)&&...
                           (poolsize>0)&&(~want_gpu);

        % Proceed with the propagation
        if parfor_makes_sense

            % Parallel loop over rho stack
            parfor k=1:numel(rho)

                % Call the commutator series procedure
                rho{k}=comm_series(spin_system,L,rho{k},time_step,nsteps);

            end

        else

            % Serial loop over rho stack
            for k=1:numel(rho)

                % Call the commutator series procedure
                rho{k}=comm_series(spin_system,L,rho{k},time_step,nsteps);

            end

        end

        % Strip the cell array if needed
        if return_numeric, rho=rho{1}; end

    end

% expm(A)*rho
else

    % Make sure state is full
    if issparse(rho), rho=full(rho); end

    % Get the scaling factor
    scaling=max(abs(rho),[],'all');

    % Catch zeros
    if scaling==0
        report(spin_system,'WARNING: all-zero state received');
        return;
    end

    % Scale the state
    rho=rho/scaling;

    % Subdivide the time step
    if isnumeric(L)
        norm_mat=cheap_norm(L)*abs(time_step);
    else
        norm_mat=max(cellfun(@cheap_norm,L))*abs(time_step);
    end
    nsteps=ceil(norm_mat/2);

    % Step checks
    if nsteps>1e4

        % Catch the common mistake of supplying wildly unreasonable |L*dt|
        error('either dt is too long, or L is too big: |L*dt|>1e4, check both.');

    elseif nsteps>100

        % Warn user if too many substeps are needed
        report(spin_system,['WARNING: ' num2str(nsteps)...
                            ' substeps required, consider using evolution() here.']);
    end

    % Decide if parallelisation would be beneficial 
    parfor_makes_sense=(size(rho,1)>256)&&(size(rho,2)>32)&&...
                       (poolsize>0)&&(~want_gpu);

    % Choose the single-vector propagation backend. The default handle below
    % is the native reordered Taylor route; the heuristic selector is called
    % only when sys.expmv_backend='auto'.
    step_backend=@(state)reordered_taylor(L,state,time_step,nsteps);
    use_backend_heuristics=isnumeric(L)&&isreal(time_step)&&(size(rho,2)==1)&&...
                           isfield(spin_system.sys,'expmv_backend')&&...
                           strcmp(spin_system.sys.expmv_backend,'auto')&&...
                           (~ismember('expv',spin_system.sys.disable));
    if use_backend_heuristics
        step_backend=step_heuristics(local_step_stats(L,time_step,norm_mat),...
                                     local_step_backends(step_backend,L,time_step));
    end

    % Proceed with the propagation
    if use_backend_heuristics

        % Call the selected single-vector backend
        rho=step_backend(rho);

    elseif parfor_makes_sense

        % Parallel loop over rho stack
        parfor m=1:size(rho,2)

            % Call Taylor-times-vector procedure
            rho(:,m)=reordered_taylor(L,rho(:,m),time_step,nsteps);

        end

    else

        % Call Taylor-times-vector procedure
        rho=reordered_taylor(L,rho,time_step,nsteps);

    end

    % Scale the result back
    rho=scaling*rho;

end

end

% Commutator series procedure
function rho=comm_series(spin_system,L,rho,time_step,nsteps)

% Convergence tolerance
tol=eps('double');

% Clean up the density matrix
rho=clean_up(spin_system,rho,tol);

% Get the scaling factor
scaling=cheap_norm(rho);

% Catch zeroes
if scaling==0
    report(spin_system,'WARNING: all-zero state received'); 
    return;
end

% Scale density matrix
rho=rho/scaling;

% Loop over substeps
for n=1:nsteps

    % Start commutator series
    next_term=rho; iter=1;

    % Sum up the series
    while nnz(next_term)>0

        % Compute the next term in the commutator series
        next_term=(-1i*(time_step/nsteps)/iter)*((L*next_term)-...
                                                 (next_term*L));

        % Clean up the term
        next_term=clean_up(spin_system,next_term,tol);

        % Add and increment the counter
        rho=rho+next_term; iter=iter+1;

    end

    % Complain if the problem is badly scaled
    if iter>32, warning(['loss of accuracy in the commutator series, iter=' num2str(iter) ...
                         ', use evolution() here instead of step()']); end

    % Final clean-up for this substep
    rho=clean_up(spin_system,rho,tol);

end

% Scale the result back
rho=scaling*rho;

end

% Taylor-times-vector procedure
function rho=reordered_taylor(L,rho,t,nsteps)

% Loop over substeps
for n=1:nsteps
    
    % Start the Taylor series
    next_term=rho; k=1; tol=eps('double');
    
    % Loop until convergence
    while nnz(abs(next_term)>tol)>0
        
        if isnumeric(L)
            
            % Centre point piecewise-constant propagator
            next_term=-(1i/k)*(t/nsteps)*(L*next_term);
            
        elseif numel(L)==2
            
            % Re-usable intermediates
            rho_a=L{1}*next_term; rho_b=L{2}*next_term;
                        
            % Left edge + right edge piecewise-linear propagator
            next_term=-(1i/2)*(1/k)*(t/nsteps)*(rho_a+rho_b)+...
                       (1/6)*(1/k)*(t^2/nsteps)*(L{1}*rho_b-L{2}*rho_a);
            
        elseif numel(L)==3
            
            % Re-usable intermediates
            rho_a=L{1}*next_term; rho_b=L{2}*next_term; rho_c=L{3}*next_term;
            
            % Left edge + midpoint + right edge piecewise-quadratic propagator
            next_term=-(1i/6)*(1/k)*(t/nsteps)*(rho_a+4*rho_b+rho_c)+...
                       (1/12)*(1/k)*(t^2/nsteps)*(L{1}*rho_c-L{3}*rho_a);
            
        end
        
        % Add the next term and increment
        rho=rho+next_term; k=k+1;
        
    end
    
    % Display non-monotonic convergence warning
    if k>32, warning(['loss of accuracy in the Taylor series, k=' num2str(k) ...
                      ', use evolution() here instead of step()']); end
    
end

end

% Step backend statistics
function stats=local_step_stats(L,time_step,norm_mat)

% These are the only data used by the heuristic selector. Keep the
% propagation state out of the policy function.
stats.matrix=L;
stats.time_step=time_step;
stats.dimension=size(L,1);
stats.is_sparse=issparse(L);
stats.is_gpu=isa(L,'gpuArray');
stats.norm_mat=norm_mat;

end

% Step backend handles
function backends=local_step_backends(default_backend,L,time_step)

% Fold the Spinach sign convention into the exponential-action backends.
% The heuristic module receives only callable handles, keeping policy
% decisions outside step() while preserving access to local subfunctions.
A=-1i*L;

% Available single-vector propagation backends. The default handle is the
% original step() Taylor path, used whenever heuristics are disabled or the
% selector chooses the legacy-like route.
backends.default=default_backend;
backends.vik=@(state)step_expmv_vik(A,state,time_step);
backends.tay1=@(state)step_expmv_tay1(A,state,time_step);
backends.tay2=@(state)step_expmv_tay2(A,state,time_step);
backends.expmv=@(state)expmv(A,state,time_step);

end

% Spinach reordered Taylor expmv backend
function w=step_expmv_vik(A,v,t)

% Estimate the generator norm
if issparse(A)&&(~isa(A,'gpuArray'))
    norm_a=normest(A,1e-2);
else
    norm_a=norm(A,1);
end

% Choose the number of substeps
nsteps=max([1 ceil(norm_a*abs(t)/2)]);

% Refuse badly scaled problems
if nsteps>1e4
    error('either dt is too long, or L is too big: |L*dt|>1e4, check both.');
end

% Get the scaling factor
scaling=max(abs(v),[],'all');

% Catch zeros
if scaling==0
    w=v; return;
end

% Scale the vector
w=v/scaling;

% Set the substep
dt=t/nsteps;

% Loop over substeps
for n=1:nsteps

    % Start the Taylor series
    term=w; k=1; tol=eps('double');

    % Loop until convergence
    while nnz(abs(term)>tol)>0

        % Taylor recurrence for the action of exp(dt*A)
        term=(dt/k)*(A*term);

        % Add the next term and increment
        w=w+term; k=k+1;

    end

end

% Scale the result back
w=scaling*w;

end

% Ibanez Algorithm 1 expmv backend
function w=step_expmv_tay1(A,v,t)

% Set algorithm parameters
m_min=40; m_max=60; s_max=45; u=2^(-53);
fact=cumprod([1 1:double(m_max+3)]);

% Fold time into the generator
if t~=1, A=t*A; end

% Get problem dimension
n=size(A,1);

% Precompute powers of the generator acting on the vector
V=zeros(n,m_max+2,'like',v);
V(:,1)=A*v;
for k=2:(m_min+2)
    V(:,k)=A*V(:,k-1);
end
n_stored=m_min+2;

% Search for Taylor degree and scaling
m=m_min; s=1; found=false();
while (~found)&&(m<=m_max)

    % Estimate the required scaling
    norm_v=norm(V(:,m+1));
    s=ceil((norm_v/(fact(m+2)*u))^(1/(m+1)));
    s=max([1 min([s_max s])]);

    % Check the two-term backward error condition
    err=norm(V(:,m+1)*(1/(s^(m+1)*fact(m+2)))+...
             V(:,m+2)*(1/(s^(m+2)*fact(m+3))));
    if err<=u
        found=true();
    else
        m=m+1;
        if m+2>n_stored
            V(:,m+2)=A*V(:,m+1);
            n_stored=m+2;
        end
    end

end

% Use the maximum degree if no pair passed the test
if ~found
    m=m_max;
    norm_v=norm(V(:,m+1));
    s=ceil((norm_v/(fact(m+2)*u))^(1/(m+1)));
    s=max([1 min([s_max s])]);
end

% Evaluate the first scaled Taylor polynomial
w=v; sk=1;
for k=1:m
    sk=sk*s;
    w=w+V(:,k)*(1/(sk*fact(k+1)));
end

% Apply the remaining scaled Taylor passes
A_scaled=A*(1/s);
for n=2:s
    z=w;
    for k=1:m
        z=A_scaled*z;
        w=w+z*(1/fact(k+1));
    end
end

end

% Ibanez Algorithm 2 expmv backend
function w=step_expmv_tay2(A,v,t)

% Set algorithm parameters
m_min=40; m_max=60; s_max=45; u=2^(-53);
fact=cumprod([1 1:double(m_max+2)]);

% Fold time into the generator
if t~=1, A=t*A; end

% Get problem dimension
n=size(A,1);

% Precompute powers of the generator acting on the vector
V=zeros(n,m_max+1,'like',v);
V(:,1)=A*v;
for k=2:(m_min+1)
    V(:,k)=A*V(:,k-1);
end
n_stored=m_min+1;

% Initialise degree and scaling search
m=m_min;
norm_v=norm(V(:,m+1));
s=max([1 min([s_max ceil((norm_v/(fact(m+2)*u))^(1/(m+1)))])]);
p=m*s; found=false();

% Search for the lowest matrix-vector product count
while (~found)&&(m<m_max)
    m_next=m+1;
    if m_next+1>n_stored
        V(:,m_next+1)=A*V(:,m_next);
        n_stored=m_next+1;
    end
    norm_next=norm(V(:,m_next+1));
    s_next=max([1 min([s_max ceil((norm_next/(fact(m_next+2)*u))^(1/(m_next+1)))])]);
    p_next=m_next*s_next;
    if p_next<=p
        m=m_next; s=s_next; p=p_next;
    else
        found=true();
    end
end

% Evaluate the first scaled Taylor polynomial
w=v; sk=1;
for k=1:m
    sk=sk*s;
    w=w+V(:,k)*(1/(sk*fact(k+1)));
end

% Apply the remaining scaled Taylor passes
A_scaled=A*(1/s);
for n=2:s
    z=w;
    for k=1:m
        z=A_scaled*z;
        w=w+z*(1/fact(k+1));
    end
end

end

% Consistency enforcement
function grumble(L,rho,time_step)
if (~isnumeric(time_step))||(~isscalar(time_step))
    error('time_step must be a scalar.');
end
if (~isnumeric(L))&&(~iscell(L))
    error('L must be a matrix or a cell array of matrices.');
end
if (~isnumeric(rho))&&(~iscell(rho))
    error('L must be a matrix, a vector, or a cell array thereof.');
end
if isnumeric(L)
    if ~allfinite(L)
        error('evolution generator is not finite.');
    end
end
if iscell(L)
    for n=1:numel(L)
        if ~allfinite(L{n})
            error('evolution generator is not finite.');
        end
    end
end
if isnumeric(rho)
    if ~allfinite(rho)
        error('state descriptor is not finite.');
    end
end
if iscell(rho)
    for n=1:numel(rho)
        if ~allfinite(rho{n})
            error('state descriptor is not finite.');
        end
    end
end
end

% Evans boldly put 50 atm of ethylene in a cell with 25 atm of oxygen. The
% apparatus subsequently blew up, but luckily not before he had obtained
% the spectra shown in Figure 8.
%
% A.J. Mehrer and R.S. Mulliken, Chem. Rev. 69 (1969) 639-656.
