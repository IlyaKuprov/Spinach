% Cayley-Magnus propagation step function. Computes the action of a
% time-dependent matrix generator on state vectors without evaluating
% matrix exponentials. Syntax:
%
%                   rho=step_cayley(spin_system,L,rho,time_step)
%
% Parameters:
%
%      L          -  Liouvillian or Hamiltonian to be used for
%                    propagation; centre point piecewise-constant
%                    rule if one matrix is supplied, fourth-order
%                    piecewise-linear rule if two matrices {left,
%                    right} are supplied, and fourth-order piecewise-
%                    quadratic rule if three matrices {left, midpoint,
%                    right} are supplied.
%
%      rho        -  state vector, or a stack of state vectors
%
%      time_step  -  length of the time step to take
%
% Outputs:
%
%      rho        -  state vector, or a stack of state vectors
%
% Notes: this function implements the five-stage fourth-order modified
%        Cayley-Magnus integrator from Blanes, Casas, and Iserles,
%        arXiv:2606.19614. Each Cayley transform is applied by BiCGSTAB,
%        with GMRES fallback if convergence fails.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=step_cayley.m>

function rho=step_cayley(spin_system,L,rho,time_step)

% Check consistency
grumble(L,rho,time_step);

% Zero time step shortcut
if time_step==0, return; end

% Do we want to run on GPU?
want_gpu=ismember('gpu',spin_system.sys.enable);

% GMRES is a CPU sparse-matrix solver
if want_gpu
    error('step_cayley.m does not support GPU execution.');
end

% Only vector Liouville/wavefunction propagation is supported here
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-wavef'})
    error('step_cayley.m only supports vector propagation formalisms.');
end

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

% Cay5_4 coefficients from Blanes-Casas-Iserles, Table 1
w31=1/(4-4^(1/3)); w21=w31; w11=1-4*w21;
w32=7/(240*(1-2*w21));
w22=(1-12*(1-w21)*w32)/(12*(1-3*w21));

% Symmetric five-stage composition
weights=[ w31  w32
          w21  w22
          w11  0
          w21 -w22
          w31 -w32 ];

% Set the linear solve tolerance
solve_tol=1e-12;

% Estimate the interval norm
if isnumeric(L)
    norm_mat=cheap_norm(L)*abs(time_step);
else
    norm_mat=max(cellfun(@cheap_norm,L))*abs(time_step);
end

% Estimate the Cayley defect tolerance
err_tol=sqrt(1e-10)/2;

% Use Spinach propagation tolerance if available
if isfield(spin_system,'tols')&&isfield(spin_system.tols,'prop_chop')&&...
   (spin_system.tols.prop_chop>0)
    err_tol=sqrt(max(spin_system.tols.prop_chop,eps('double')))/2;
end

% Estimate the leading fourth-order Cayley rational defect
if isnumeric(L)
    err_coeff=abs(sum(weights(:,1).^5))/80;

    % Apply the leading defect operator to the state stack
    err_state=rho;
    for k=1:5
        err_state=(-1i*time_step)*(L*err_state);
    end

    % Estimate the relative Cayley defect on the supplied state
    state_norm=max(norm(rho,'fro'),eps('double'));
    cayley_err=err_coeff*norm(err_state,'fro')/state_norm;
else
    wgt_norm=abs(weights(:,1))+2*abs(weights(:,2));
    err_coeff=sum(wgt_norm.^5)/80;

    % Use a norm bound for sampled time-dependent generators
    cayley_err=err_coeff*norm_mat^5;
end

% Estimate the Cayley stage norm
if isnumeric(L)
    stage_norm=max(abs(weights(:,1)))*norm_mat/2;
else
    stage_norm=max(wgt_norm)*norm_mat/2;
end

% Choose the subdivision count from accuracy and conditioning
nsteps=max([1 ceil((cayley_err/err_tol)^(1/4)) ceil(stage_norm)]);

% Step checks
if norm_mat>1e4

    % Catch the common mistake of supplying wildly unreasonable |L*dt|
    error('either dt is too long, or L is too big: |L*dt|>1e4, check both.');

elseif nsteps>100

    % Warn user if too many substeps are needed
    report(spin_system,['WARNING: ' num2str(nsteps)...
                        ' substeps required, consider using evolution() here.']);

elseif iscell(L)&&(norm_mat>100)

    % Warn user if sampled interval is too long
    report(spin_system,'WARNING: large sampled interval, consider shorter external time steps.');

end

% Propagate the vector stack
rho=cayley_magnus(L,rho,time_step,nsteps,weights,solve_tol);

% Scale the result back
rho=scaling*rho;

end

% Fourth-order modified Cayley-Magnus integrator
function rho=cayley_magnus(L,rho,t,nsteps,weights,solve_tol)

% Loop over substeps
for n=1:nsteps

    % Resample the supplied interpolation model
    if iscell(L)&&(nsteps>1)

        % Get the local left interval position
        pos_l=(n-1)/nsteps;

        % Get the local right interval position
        pos_r=n/nsteps;

        % Piecewise-linear interpolation
        if numel(L)==2

            % Interpolate the local endpoint generators
            L_sub={L{1}+pos_l*(L{2}-L{1}),L{1}+pos_r*(L{2}-L{1})};

        % Piecewise-quadratic interpolation
        elseif numel(L)==3

            % Get the local midpoint position
            pos_m=(n-1/2)/nsteps;

            % Interpolate the local three-point generators
            L_sub={((2*pos_l^2-3*pos_l+1)*L{1})+...
                   ((4*pos_l-4*pos_l^2)*L{2})+...
                   ((2*pos_l^2-pos_l)*L{3}),...
                   ((2*pos_m^2-3*pos_m+1)*L{1})+...
                   ((4*pos_m-4*pos_m^2)*L{2})+...
                   ((2*pos_m^2-pos_m)*L{3}),...
                   ((2*pos_r^2-3*pos_r+1)*L{1})+...
                   ((4*pos_r-4*pos_r^2)*L{2})+...
                   ((2*pos_r^2-pos_r)*L{3})};

        end

    else

        % Use the supplied generator representation directly
        L_sub=L;

    end

    % Generate Cayley-Magnus basis matrices
    [alpha1,alpha2]=alphas(L_sub,t/nsteps);

    % Apply the Cayley maps
    for k=size(weights,1):-1:1
        rho=apply_cayley(weights(k,1)*alpha1+weights(k,2)*alpha2,rho,solve_tol);
    end

end

end

% Basis matrices for the fourth-order scheme
function [alpha1,alpha2]=alphas(L,t)

% Spinach uses exp(-1i*L*t), paper uses exp(A*t)
if isnumeric(L)

    % Centre point piecewise-constant generator
    alpha1=-1i*t*L; alpha2=0*L;

elseif numel(L)==2

    % Fourth-order piecewise-linear rule
    alpha1=(-1i*t/2)*(L{1}+L{2});
    alpha2=(-1i*t)*(L{2}-L{1});

elseif numel(L)==3

    % Fourth-order piecewise-quadratic rule
    alpha1=(-1i*t/6)*(L{1}+4*L{2}+L{3});
    alpha2=(-1i*t)*(L{3}-L{1});

else

    % Unsupported cell shape
    error('L must be a matrix, a two-matrix cell array, or a three-matrix cell array.');

end

end

% Cayley transform action, (I-X/2)\(I+X/2)*rho
function rho=apply_cayley(X,rho,tol)

% Common parameters
maxit=min(size(X,1),max(32,ceil(size(X,1)/4)));
A=speye(size(X,1))-X/2; B=speye(size(X,1))+X/2;

% Precompute the right-hand side stack
rhs=B*rho; next_rho=zeros(size(rhs),'like',rhs);

% Loop over right-hand sides
for n=1:size(rho,2)

    % BiCGSTAB solve
    [next_rho(:,n),flag]=bicgstab(A,rhs(:,n),tol,maxit,[],[],rhs(:,n));

    % Fall back to GMRES on solver failure
    if flag~=0
        [next_rho(:,n),flag]=gmres(A,rhs(:,n),[],tol,maxit,[],[],rhs(:,n));
    end

    % Report solver failure
    if flag~=0
        error('linear solver failed to converge in a Cayley-Magnus stage.');
    end

end

% Return the solved stack
rho=next_rho;

end

% Consistency enforcement
function grumble(L,rho,time_step)
if (~isnumeric(time_step))||(~isscalar(time_step))||(~isfinite(time_step))
    error('time_step must be a finite scalar.');
end
if (~isnumeric(L))&&(~iscell(L))
    error('L must be a matrix or a cell array of matrices.');
end
if ~isnumeric(rho)
    error('rho must be a vector or a matrix of column-stacked vectors.');
end
if isnumeric(L)
    if ~allfinite(L)
        error('evolution generator is not finite.');
    end
    if size(L,1)~=size(L,2)
        error('evolution generator must be square.');
    end
end
if iscell(L)
    if ~ismember(numel(L),[2 3])
        error('L cell array must have two or three matrices.');
    end
    for n=1:numel(L)
        if ~isnumeric(L{n})||(~allfinite(L{n}))
            error('evolution generator is not finite.');
        end
        if size(L{n},1)~=size(L{n},2)
            error('evolution generator must be square.');
        end
    end
end
if ~allfinite(rho)
    error('state descriptor is not finite.');
end
end
