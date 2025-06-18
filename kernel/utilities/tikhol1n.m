% L1 norm Tikhonov regularised solver for A*x=y where
% A is an ill-conditioned matrix. The error functional
% is norm(A*x-y,1)^2+lambda*norm(x,1), it is minimised
% using the FISTA algorithm. Syntax:
%
%           [x,err,reg]=tikhol1n(A,y,lambda)
%
% Parameters:
%
%    A - a real or complex matrix
%
%    y - a real or complex column vector
%
%    lambda - a non-negative real scalar
%
% Outputs:
%
%    x - a real or complex vector
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=tikhol1n.m>

function [x,err,reg]=tikhol1n(A,y,lambda)

% Tolerances
normest_tol=1e-3;   % relative 2-norm estimation tolerance
step_norm_tol=1e-5; % relative step norm convergence tolerance
rel_nz_tol=1e-6;    % relative (to the max) non-zero tolerance

% Pre-compute CT
A_ct=ctranspose(A);
        
% Lipschitz constant and the threshold
L=2*normest(A,normest_tol)^2; thr=lambda/L;

% Complex soft thresholding function
soft_thr=@(x)sign(x).*max(abs(x)-thr,0);
        
% Complex random initial guess
x=   randn(size(A,2),1)+...
  1i*randn(size(A,2),1); x_old=x;

% Iteration counters and momentum
t=1; iter_count=0; converged=false;

% Report to the user
disp(['FISTA called with lambda = ' num2str(lambda)]);
        
% FISTA iteration loop
while ~converged

    % Get the error vector
    err_vec=A*x-y;

    % Get the gradient of the error
    g=2*(A_ct*err_vec);

    % Compute the proximal point
    x_prox=soft_thr(x-g/L);

    % Update Nesterov momentum
    t_new=0.5*(1+sqrt(1+4*t^2));

    % Take the step
    x=x_prox+((t-1)/t_new)*(x_prox-x_old);

    % Check convergence
    soln_norm=norm(x_prox,2);
    step_norm=norm(x_prox-x_old,2);
    converged=(step_norm<step_norm_tol*soln_norm);

    % Close the loop
    t=t_new; x_old=x_prox;
    iter_count=iter_count+1;

    % Compute figures of merit
    err=norm(err_vec,2)^2; reg=norm(x,1);

    % Report progress
    if mod(iter_count,100)==0
        zf=nnz(abs(x)<rel_nz_tol*max(abs(x)))/numel(x);
        disp(['FISTA iter ' int2str(iter_count) ', zf ' num2str(zf) ...
              ', err '  num2str(err) ', reg ' num2str(reg)          ...
              ', step ' num2str(step_norm/soln_norm)]);
    end

end

end

% "He begins working on calculus problems in his head as 
% soon as he awakens. He did calculus while driving in his
% car, while sitting in the living room, and while lying 
% in bed at night."
%
% Mary Louise Bell (Richard Feynman's second
% wife), in her divorce complaint

