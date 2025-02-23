% Propagation step function. Uses Taylor propagation and sparse exponenti-
% ation where appropriate. Supports piecewise-constant, piecewise-linear,
% and piecewise quadratic product quadratures. Syntax:
%
%                   rho=step(spin_system,L,rho,time_step)
%
% Arguments:
%
%      L          -  Liouvillian or Hamiltonian to be used for 
%                    propagation; centre point piecewise-constant
%                    rule if one matrix is supplied, piecewise-
%                    linear rule if two matrices {left, right} 
%                    are supplied, piecewise-quadratic if three
%                    matrices {left, midpoint, right} are given.
%
%      rho        -  state vector or density matrix to be propagated
%
%      time_step  -  length of the time step to take
%
% Note: we initially had a faithful implementation of the Krylov process
%       here - subspace, orthogonalisation, projection, etc., but in all
%       our testing it was much inferior to the reordered Taylor process
%       that is currently implemented below.
%
% Note: the peculiar sequence of algebraic operations in the code below
%       is designed to minimise the memory footprint in large cases.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
% a.acharya@soton.ac.uk
% c.musselwhite@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=step.m>

function rho=step(spin_system,L,rho,time_step)

% Check consistency
grumble(L,rho,time_step);

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

    % Fast serial bypass for small matrices
    if size(L,1)<spin_system.tols.small_matrix

        % Use Matlab's expm
        P=expm(-1i*L*time_step);
        rho=P*rho*P'; return;

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

    % Proceed with the propagation
    if parfor_makes_sense

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

% Consistency enforcement
function grumble(L,rho,time_step)
if (~isnumeric(time_step))||(~isscalar(time_step))
    error('time_step must be a scalar.');
end
if (~isnumeric(L))&&(~iscell(L))
    error('L must be a matrix or a cell array of matrices');
end
if (~isnumeric(rho))&&(~iscell(rho))
    error('L must be a matrix, a vector, or a cell array thereof');
end
end

% Evans boldly put 50 atm of ethylene in a cell with 25 atm of oxygen. The
% apparatus subsequently blew up, but luckily not before he had obtained
% the spectra shown in Figure 8.
%
% A.J. Mehrer and R.S. Mulliken, Chem. Rev. 69 (1969) 639-656.

