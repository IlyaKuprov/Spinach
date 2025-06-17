% L1 norm Tikhonov regularised solver for ...
%
% ilya.kuprov@weizmann.ac.il

function [spec,err,reg]=tikhol1n(A,y,lambda)

% Tolerances
normest_tol=1e-3;
step_norm_tol=1e-4;

% Pre-compute CT
A_ct=ctranspose(A);
        
% Lipschitz constant and the threshold
L=2*normest(A,normest_tol)^2; thr=lambda/L;

% Complex soft thresholding function
soft_thr=@(x)sign(x).*max(abs(x)-thr,0);
        
% Complex random initial guess
spec=   randn(size(A,2),1)+...
     1i*randn(size(A,2),1);
spec_old=spec;

% Iteration counters and momentum
t=1; iter_count=0; converged=false;
        
% FISTA iteration loop
while ~converged

    % Get the error vector
    err_vec=A*spec-y;

    % Get the gradient of the error
    g=2*(A_ct*err_vec);

    % Compute the proximal point
    spec_prox=soft_thr(spec-g/L);

    % Update Nesterov momentum
    t_new=0.5*(1+sqrt(1+4*t^2));

    % Take the step
    spec=spec_prox+((t-1)/t_new)*(spec_prox-spec_old);

    % Check convergence
    soln_norm=norm(spec_prox,2);
    step_norm=norm(spec_prox-spec_old,2);
    converged=(step_norm<step_norm_tol*soln_norm);

    % Close the loop
    t=t_new; spec_old=spec_prox;
    iter_count=iter_count+1;

    % Compute figures of merit
    err=norm(err_vec,2)^2; reg=norm(spec,1);

    % Report progress
    if mod(iter_count,10)==0
        disp(['FISTA iter ' int2str(iter_count) ...
              ', err '  num2str(err) ', reg ' num2str(reg) ...
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

