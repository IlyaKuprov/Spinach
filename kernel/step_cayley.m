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

% Estimate the step norm
if isnumeric(L)
    norm_mat=cheap_norm(L)*abs(time_step);
else
    norm_mat=max(cellfun(@cheap_norm,L))*abs(time_step);
end

% Sampled time-dependent intervals cannot be subdivided without resampling
if isnumeric(L)
    nsteps=ceil(norm_mat/2);
else
    nsteps=1;
end

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
rho=cayley_magnus(L,rho,time_step,nsteps);

% Scale the result back
rho=scaling*rho;

end

% Fourth-order modified Cayley-Magnus integrator
function rho=cayley_magnus(L,rho,t,nsteps)

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

% Loop over substeps
for n=1:nsteps

    % Generate Cayley-Magnus basis matrices
    [alpha1,alpha2]=alphas(L,t/nsteps);

    % Apply the Cayley maps
    for k=size(weights,1):-1:1
        rho=apply_cayley(weights(k,1)*alpha1+weights(k,2)*alpha2,rho);
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
function rho=apply_cayley(X,rho)

% Common parameters
tol=1e-12; maxit=min(size(X,1),max(32,ceil(size(X,1)/4)));
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
